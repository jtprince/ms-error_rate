require 'ms/fasta'
require 'ms/in_silico/digester'

module Ms
  module Ident

    IPI_RE = /IPI:([\w\d\.]+)\|/
    GI_RE = /gi|([\w\d\.]+)\|/

    # the twenty standard amino acids
    STANDARD_AA = %w(A C D E F G H I K L M N P Q R S T V W Y)

    DEFAULT_PEPTIDE_CENTRIC_DB = {:missed_cleavages => 1, :min_length => 8, :enzyme => Ms::InSilico::Digester::TRYPSIN, :id_regexp => nil, :remove_digestion_file => true, :cleave_initiator_methionine => true, :expand_aa => {'X' => STANDARD_AA}}

    # writes a new file with the added 'min_aaseq<Integer>'
    # creates a temporary digestion file that contains all peptides digesting
    # with certain missed_cleavages (i.e., min_seq_length is not applied to
    # this file but on the final peptide centric db)
    def self.peptide_centric_db(fasta_file, opts={})
      opts = DEFAULT_PEPTIDE_CENTRIC_DB.merge(opts)

      (missed_cleavages, min_length, enzyme, id_regexp, remove_digestion_file, cleave_initiator_methionine, expand_aa) = opts.values_at(:missed_cleavages, :min_length, :enzyme, :id_regexp, :remove_digestion_file, :cleave_initiator_methionine, :expand_aa) 

      unless id_regexp
        id_regexp = Ms::Fasta.id_regexp(Ms::Fasta.filetype(fasta_file))
        raise RuntimeError, "fasta file type not recognized, supply id_regexp" unless id_regexp
      end

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
        id = prot.match(id_regexp)[1]
        peps.each do |pep|
          if pep.size >= min_length
            hash[pep] << id
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
          out.puts( "#{k}: #{v.join('-')}" )
        end
      end
      puts "#{Time.now - start_time} sec" if $VERBOSE

      if remove_digestion_file
        File.unlink(digestion_file)
      end

    end

    # does combinatorial expansion of all letters requesting it.
    # expand_aa is hash like: {'X'=>STANDARD_AA}
    def self.expand_peptides(peptide, expand_aa)
      letters_in_order = expand_aa.keys.sort
      index_and_key = []
      peptide.split('').each_with_index do |char,i| 
        if let_index = letters_in_order.index(char) 
          index_and_key << [i, letters_in_order[let_index]]
        end
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
