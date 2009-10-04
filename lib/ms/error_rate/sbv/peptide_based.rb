require 'ms/error_rate/sbv'

module Ms
  module ErrorRate
    class Sbv
      # Constraints on aaseq attribute of peptides (the bare amino acid sequence)
      # works by calculating amino acid frequencies in the fasta file used.
      class PeptideBased

        def self.generate_hashes(pep_to_prot_file, aa="C", min_num=1 )
          Ms::ErrorRate::Sbv.generate_hashes(pep_to_prot_file, :type_code => "#{TWO_LETTER}_min#{min_num}") do |pep|
            if min_num == 1
              if pep.include?(aa) ; 1
              else ; 0
              end
            else
              count = 0
              pep.each_char {|c| count += 1 if c == aa }
              if count >= min_num ; 1
              else ; 0
              end
            end
          end
        end

      end # class
    end # Sbv
  end # ER
end # Ms

