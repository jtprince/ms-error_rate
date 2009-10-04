

module Ms
  module ErrorRate
    # Sample Bias Validator
    class Sbv
      LENGTH_EXT = 'freq_by_length'
      AASEQ_EXT = 'by_aaseq'

      # if a protein hash is given, will yield the return an array of
      # values generated with the value from keying each protein of the
      # peptide.  Otherwise, will yield each peptide in turn
      def self.generate_hashes(pep_to_prot_file, opts={})
        op = { :aaseq_ext => AASEQ_EXT,
          :length_ext => LENGTH_EXT,
          :file_ext => '.yml',
          :type_code => '',
          :protein_hash => nil, 
          :stderr_counter => true,
        }.merge(opts)

        base = pep_to_prot_file.chomp(File.extname(pep_to_prot_file)) 
        freqs = Hash.new {|h,k| h[k] = 0.0 }
        counts = Hash.new {|h,k| h[k] = 0 }
        (fileout1, fileout2) = [:aaseq_ext, :length_ext].map do |type_ext|
          base + '.' + op[:type_code] + '.' + op[type_ext] + op[:file_ext]
        end
        protein_hash = op[:protein_hash]
        pep_count = 0
        if op[:stderr_counter]
          $stderr.print "[working, 100,000 peptides = '.'] "
          $stderr.flush
        end
        File.open(fileout1 , 'w') do |out|
          IO.foreach(pep_to_prot_file) do |line|
            (pep, prot_string) = line.chomp!.split(': ')

            total_transmembrane = 0 
            total_known = 0
            answ = 
              if protein_hash
                yield( protein_hash.values_at(prot_string.split('-')) )
              else
                yield(pep)
              end
            out.puts "#{pep}: #{answ}"
            freqs[pep.size] += answ
            counts[pep.size] += 1
            pep_count += 1
            if pep_count % 100000 == 0 && op[:stderr_counter]
              $stderr.print '.'
              $stderr.flush
            end
          end
        end
        $stderr.puts "DONE!" if op[:stderr_counter]
        avg_freq_ar = {}
        freqs.each do |k,v|
          avg_freq_ar[k] = v / counts[k]
        end
        File.open(fileout2, 'w') {|out| out.print avg_freq_ar.to_yaml }
        [fileout1, fileout2] 
      end


      # a hash by aaseq giving a value between 0 and 1 telling how much of
      # an indicator the hit is
      attr_accessor :indicator_by_aaseq

      attr_accessor :frequency_indicator_opposite

      attr_accessor :size_to_freq

      # boolean
      attr_accessor :indicators_signify_true_hit


      def initialize(indicator_by_aaseq_hash, size_to_freq, frequency_indicator_opposite, indicators_signify_true_hit=true)
        @indicators_signify_true_hit = indicators_signify_true_hit
        @frequency_indicator_opposite = frequency_indicator_opposite
        @indicator_by_aaseq = indicator_by_aaseq_hash
        @tot_num_indicators = 0.0
        @tot_num = 0
      end

      # returns the cumulative precision (fraction of true positives among
      # total hits) frequency_of_indicators is the probability that a generic
      # amino acid sequence will be an indicator (this may variable by
      # sequence length).
      def update_precision(aaseq)
        @tot_num_indicators << indicator_by_aaseq[aaseq]
        @tot_num += 1
        @frequency_of_indicators_sum += @size_to_freq[aaseq.size]
        # FP Indicator
        value = @tot_num_indicators * (1.0 - @frequency_indicator_opposite) * @frequency_of_indicators_sum / (@tot_num**2)
        precision = 
          if @indicators_signify_true_hit 
            value  # a true indicator type (gives precision)
          else  # false indicator type
            1 - value  # 1 - fdr == precision
          end
      end

      def calculate_background_frequency
        @aaseq_to_fraction
      end

    end

  end
end
