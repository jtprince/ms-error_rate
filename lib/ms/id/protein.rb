require 'trie'

module Ms ; end
module Ms::Id ; end

module Ms::Id::Protein

  # takes a group of peptide hits and a simple_protein_digested_file
  # each line has (protein seq seq seq ...)
  def self.minimal_protein_set(hits, prots_digested_file)

    # this is a very fast trie but cannot hold arrays, just strings or
    # numbers!
    tr = Trie.new 
    prot_num = 0
    IO.foreach(prots_digested_file) do |line|
      (prot, *peps) = line.chomp!.split(/\s+/)
      peps.each do |pep|
        if pep.size >= 4
          if tr.has_key?(pep)
            p pep
            puts "HERE1"
            puts (tr.get(pep) + ",#{prot_num}")
            tr.add pep, (tr.get(pep) + ",#{prot_num}")
          else
            puts "HERE2"
            tr.add pep, prot_num.to_s
          end
        end
      end
      prot_num += 1
    end
    tr
  end

end
