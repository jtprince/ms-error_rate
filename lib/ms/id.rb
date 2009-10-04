require 'ms/fasta'
require 'ms/in_silico/digester'

module Ms
  module Id

    IPI_RE = /IPI:([\w\d\.]+)\|/
    GI_RE = /gi|([\w\d\.]+)\|/

    DEFAULT_PEPTIDE_CENTRIC_DB = {:missed_cleavages => 1, :min_length => 8, :enzyme => Ms::InSilico::Digester::TRYPSIN, :id_regexp => nil, :remove_digestion_file => true }

    # writes a new file with the added 'min_aaseq<Integer>'
    # creates a temporary digestion file that contains all peptides digesting
    # with certain missed_cleavages (i.e., min_seq_length is not applied to
    # this file but on the final peptide centric db)
    def self.peptide_centric_db(fasta_file, opts={})
      opts = DEFAULT_PEPTIDE_CENTRIC_DB.merge(opts)

      (missed_cleavages, min_length, enzyme, id_regexp, remove_digestion_file) = opts.values_at(:missed_cleavages, :min_length, :enzyme, :id_regexp, :remove_digestion_file) 

      unless id_regexp
        id_regexp = Ms::Fasta.id_regexp(Ms::Fasta.filetype(fasta_file))
        raise RuntimeError, "fasta file type not recognized, supply id_regexp" unless id_regexp
      end

      base = fasta_file.chomp(File.extname(fasta_file))
      digestion_file = base + ".msd_clvg#{missed_cleavages}.peptides"
      File.open(digestion_file, "w") do |fh|
        peptides = []
        Ms::Fasta.open(fasta_file) do |fasta|
          fasta.each do |prot|
            fh.puts( prot.header.split(/\s+/).first + "\t" + enzyme.digest(prot.sequence, missed_cleavages).join(" ") )
          end
        end
      end

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
      base = digestion_file.chomp(File.extname(digestion_file))
      File.open(base + ".min_aaseq#{min_length}" + ".yml", 'w') do |out|
        hash.each do |k,v|
          out.puts( "#{k}: #{v.join('-')}" )
        end
      end

      if remove_digestion_file
        File.unlink(digestion_file)
      end

    end

  end
end
