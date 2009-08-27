

module Ms
  module ErrorRate
    module SampleBias
      # peptides must respond to transmembrane and return a value of 0 to 1.
      # 0 means all its associated proteins do not have tm segments.  1 means
      # all of them do.  Fractions OK.
      def transmembrane(peptides, 
    end
  end
end
