require 'ms/fasta'
require 'ms/in_silico/digester'

module Ms
  module Id

    IPI_RE = /IPI:([\w\d\.]+)\|/
    GI_RE = /gi|([\w\d\.]+)\|/

    DEFAULT_PEPTIDE_CENTRIC_DB = {:missed_cleavages => 1, :min_length => 8, :enzyme => Ms::InSilico::Digester::TRYPSIN, :id_regexp => nil, :remove_digestion_file => true, :also_cleave_initiator_methionine => true}

    # writes a new file with the added 'min_aaseq<Integer>'
    # creates a temporary digestion file that contains all peptides digesting
    # with certain missed_cleavages (i.e., min_seq_length is not applied to
    # this file but on the final peptide centric db)
    def self.peptide_centric_db(fasta_file, opts={})
      opts = DEFAULT_PEPTIDE_CENTRIC_DB.merge(opts)

      (missed_cleavages, min_length, enzyme, id_regexp, remove_digestion_file, also_cleave_initiator_methionine) = opts.values_at(:missed_cleavages, :min_length, :enzyme, :id_regexp, :remove_digestion_file, :also_cleave_initiator_methionine) 

      unless id_regexp
        id_regexp = Ms::Fasta.id_regexp(Ms::Fasta.filetype(fasta_file))
        raise RuntimeError, "fasta file type not recognized, supply id_regexp" unless id_regexp
      end

      start_time = Time.now
      print "Digesting #{fasta_file} ..." if $VERBOSE

      base = fasta_file.chomp(File.extname(fasta_file))
      digestion_file = base + ".msd_clvg#{missed_cleavages}.peptides"
      File.open(digestion_file, "w") do |fh|
        Ms::Fasta.open(fasta_file) do |fasta|
          fasta.each do |prot|
            peptides = enzyme.digest(prot.sequence, missed_cleavages)
            begin_peptides_TEST = peptides.map
            if (also_cleave_initiator_methionine && (prot.sequence[0,1] == "M"))
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

  end
end
