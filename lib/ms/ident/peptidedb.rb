require 'ms/fasta'
require 'ms/in_silico/digester'

module Ms
  module Ident

    class Peptidedb
      MAX_NUM_AA_EXPANSION = 3

      # the twenty standard amino acids
      STANDARD_AA = %w(A C D E F G H I K L M N P Q R S T V W Y)

      DEFAULT_PEPTIDE_CENTRIC_DB = {:missed_cleavages => 2, :min_length => 4, :enzyme => Ms::InSilico::Digester::TRYPSIN, :id_regexp => nil, :remove_digestion_file => true, :cleave_initiator_methionine => true, :expand_aa => {'X' => STANDARD_AA}}

      def self.cmdline(argv)

        opt = {
          :remove_digestion_file => true,
          :enzyme => Ms::InSilico::Digester::TRYPSIN
        }
        opts = OptionParser.new do |op|
          op.banner = "usage: #{File.basename($0)} <file>.fasta ..."
          op.separator "output: "
          op.separator "    <file>.msd_clvg<missed_cleavages>.min_aaseq<min_length>.yml"
          op.separator "format:"
          op.separator "    PEPTIDE: ID1<tab>ID2<tab>ID3..."
          op.separator ""    
          op.separator "    Initiator Methionines - by default, will generate two peptides"
          op.separator "    for any peptide found at the N-termini starting with 'M'"
          op.separator "    (i.e., one with and one without the leading methionine)"
          op.separator ""
          op.on("--missed-cleavages <#{opt[:missed_cleavages]}>", Integer, "max num of missed cleavages") {|v| opt[:missed_cleavages] = v }
          op.on("--min-length <#{opt[:min_length]}>", Integer, "the minimum peptide aaseq length") {|v| opt[:min_length] = v }
          op.on("--no-cleaved-methionine", "does not cleave off initiator methionine") { opt[:cleave_initiator_methionine] = false }
          op.on("--no-expand-x", "don't enumerate aa 'X' possibilities") { opt[:expand_aa] = nil }
          op.on("-e", "--enzyme <name>", "enzyme for digestion") {|v| opt[:enzyme] = Ms::Insilico::Digester.const_get(v.upcase) }
          op.on("--list-enzymes", "lists approved enzymes and exits") do
            puts Ms::InSilico::Digester::ENZYMES.keys.join("\n")
            exit
          end
        end

        opts.parse!(argv)

        if argv.size == 0
          puts opts || exit
        end

        argv.each do |file|
          Ms::Ident::Peptidedb.peptide_centric_db(file, opt)
        end
      end

      # writes a new file with the added 'min_aaseq<Integer>'
      # creates a temporary digestion file that contains all peptides digesting
      # with certain missed_cleavages (i.e., min_seq_length is not applied to
      # this file but on the final peptide centric db)
      def self.peptide_centric_db(fasta_file, opts={})
        opts = DEFAULT_PEPTIDE_CENTRIC_DB.merge(opts)

        (missed_cleavages, min_length, enzyme, id_regexp, remove_digestion_file, cleave_initiator_methionine, expand_aa) = opts.values_at(:missed_cleavages, :min_length, :enzyme, :id_regexp, :remove_digestion_file, :cleave_initiator_methionine, :expand_aa) 
        start_time = Time.now
        print "Digesting #{fasta_file} ..." if $VERBOSE

        if expand_aa
          letters_to_expand_re = Regexp.new("[" << Regexp.escape(expand_aa.keys.join) << "]")
        end

        base = fasta_file.chomp(File.extname(fasta_file))
        digestion_file = base + ".msd_clvg#{missed_cleavages}.peptides"
        File.open(digestion_file, "w") do |fh|
          Ms::Fasta.open(fasta_file) do |fasta|
            fasta.each do |prot|
              peptides = enzyme.digest(prot.sequence, missed_cleavages)
              if (cleave_initiator_methionine && (prot.sequence[0,1] == "M"))
                m_peps = []
                init_methionine_peps = []
                peptides.each do |pep|
                  # if the peptide is at the beginning of the protein sequence
                  if prot.sequence[0,pep.size] == pep
                    m_peps << pep[1..-1]
                  end
                end
                peptides.push(*m_peps)
              end
              if expand_aa
                peptides = peptides.map do |pep|
                  if pep =~ letters_to_expand_re
                    expand_peptides(pep, expand_aa)
                  else
                    pep
                  end
                end.flatten
              end
              fh.puts( prot.header.split(/\s+/).first + "\t" + peptides.join(" ") )
            end
          end
        end
        puts "#{Time.now - start_time} sec" if $VERBOSE


        start_time = Time.now
        print "Organizing raw digestion #{digestion_file} ..." if $VERBOSE

        hash = Hash.new {|h,k| h[k] = [] }
        IO.foreach(digestion_file) do |line|
          (prot, *peps) = line.chomp!.split(/\s+/)
          # prot is something like this: "sp|P31946|1433B_HUMAN" in uniprot 
          peps.each do |pep|
            if pep.size >= min_length
              hash[pep] << prot
            end
          end
        end
        puts "#{Time.now - start_time} sec" if $VERBOSE

        base = digestion_file.chomp(File.extname(digestion_file))
        final_outfile = base + ".min_aaseq#{min_length}" + ".yml"

        start_time = Time.now
        print "Writing results to #{} ..." if $VERBOSE

        File.open(final_outfile, 'w') do |out|
          hash.each do |k,v|
            out.puts( [k, v.join("\t")].join(": ") )
          end
        end
        puts "#{Time.now - start_time} sec" if $VERBOSE

        if remove_digestion_file
          File.unlink(digestion_file)
        end

      end

      # does combinatorial expansion of all letters requesting it.
      # expand_aa is hash like: {'X'=>STANDARD_AA}
      # returns nil if there are more than MAX_NUM_AA_EXPANSION amino acids to
      # be expanded
      # returns an empty array if there is no expansion
      def self.expand_peptides(peptide, expand_aa)
        letters_in_order = expand_aa.keys.sort
        index_and_key = []
        peptide.split('').each_with_index do |char,i| 
          if let_index = letters_in_order.index(char) 
            index_and_key << [i, letters_in_order[let_index]]
          end
        end
        if index_and_key.size > MAX_NUM_AA_EXPANSION
          return nil
        end
        to_expand = [peptide]
        index_and_key.each do |i,letter|
          new_peps = []
          while current_pep = to_expand.shift do
            new_peps << expand_aa[letter].map {|v| dp = current_pep.dup ; dp[i] = v ; dp }
          end
          to_expand = new_peps.flatten
        end
        to_expand
      end

    end
  end
end
