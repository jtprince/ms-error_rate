
module Ms
  module ErrorRate
    module Decoy
      # this is the # true positives (found by estimating the number of false
      # hits using the # decoy)
      # frit == fraction 
      def self.precision(num_target, num_decoy, frit=1.0)
        # will calculate as floats in case fractional amounts passed in for
        # whatever reason
        num_target_f = num_target.to_f
        num_true_pos = num_target_f - (num_decoy.to_f * frit)
        precision =
          if num_target_f == 0.0
            if num_decoy.to_f > 0.0
              0.0
            else
              1.0
            end
          else
            num_true_pos/num_target_f
          end
        precision
      end
    end
  end
end
