#require 'trie'

module Containers ; end

require 'rubygems'
require 'containers/trie'

module Ms ; end
module Ms::Id ; end

module Ms::Id::Protein

  # takes a group of peptide hits and a simple_protein_digested_file
  # each line has (protein seq seq seq ...)
  def self.minimal_protein_set(hits, prots_digested_file)

    # this is a very fast trie but cannot hold arrays, just strings or
    # numbers!
    #tr = Trie.new 
    tr = Containers::Trie.new
    #prot_num = 0
    IO.foreach(prots_digested_file) do |line|
      (prot, *peps) = line.chomp!.split(/\s+/)
      ipi = prot.match(/IPI:([\w\d\.]+)\|/)[1]
      peps.each do |pep|
        if pep.size >= 9

          if x = tr[pep]
            #puts (tr.get(pep) + ",#{prot_num}")
            tr[pep] =  x + ",#{ipi}"
          else
            tr[pep] = ipi
          end

          #if tr.has_key?(pep)
          #  #puts (tr.get(pep) + ",#{prot_num}")
          #  tr.add pep, (tr.get(pep) + "-#{prot_num}")
          #else
          #  tr.add pep, prot_num.to_s
          #end


        end
      end
      #prot_num += 1
    end
    tr
  end

end
