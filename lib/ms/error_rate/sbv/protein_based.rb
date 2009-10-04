require 'ms/fasta'
require 'ms/error_rate/sbv'
require 'transmembrane'

module Ms
  module ErrorRate
    class Sbv
      module ProteinBased
        # note the pep to prot hash has proteins in a string separated by a
        # hyphen.  returns the names of the files written
        def self.generate_hashes(pep_to_prot_file, protid_to_val, options={})
          options[:protein_hash] = protid_to_val
          files = Ms::ErrorRate::Sbv.generate_hashes(pep_to_prot_file, options) do |prot_return_vals|
            
            total_with_bias = 0 
            total_known = 0
            prot_return_vals.each do |val|
              if !val.nil?
                total_with_bias += val
                total_known += 1
              end
            end
            total_with_bias.to_f / total_known
          end #block

          files

        end # end method
      end # module
    end # class
  end # ErrorRate
end # Ms 

