
require 'ms/error_rate/decoy'

module Ms
  module ErrorRate
    # For generating and working with q-value calculations.  The q-value is the global false discovery rate when accepting that particular ID.  We do not necessarily distinguish here between *how* the FDR is generated (i.e., Storey's pFDR "the occurrence of false positives" vs. Benjamini-Hochberg's FDR "the rate of false positives" [except to prefer Storey when possible] ).  The main point is that we sort and threshold based on a global FDR.
    module Qvalue


         
      # returns a parallel array to target hits with qvalues
      # opts = :z_together true/false (default false) group all charges
      # together.
      # the sort block should sort from worst to best
      # by default, sorting is: {|hit| hit.score} if not provided
      def target_decoy_qvalues(target_hits, decoy_hits, opts={}, &sorting=nil)
        sorting = lambda {|hit| hit.score } unless sorting
        opts = {:z_together => false}.merge(opts)
        target_set = Set.new(target_hits)

        # Proc.new doesn't do arity checking
        qvalues = Proc.new do |hits|
          qvalues = []
          sorted_best_to_worst = (hits.sort_by &sorting).reverse
          (target_hits, qvalues) = Ms::Id::Qvalue.mixed_target_decoy(sorted_best_to_worst.reverse, target_set)
          target_hits.zip(qvalues) do |hit, qvalue|
            qvalues << qvalue
          end
          qvalues
        end

        if !opts[:z_together]
          target_hits
        end

        if opts[:z_together]
          all_hits = []
          all_hits_by_charge.each do |charge, hits|
            all_hits.push(*hits)
          end
          hits_with_qvalues.call(all_hits)
        else
          all_hits = []
          all_hits_by_charge.each do |charge,hits|
            all_hits.push(*(hits_with_qvalues.call(hits)))
          end
          all_hits.sort_by {|v| -(v.last) }
        end

      end

      # returns [target_hits, qvalues] (parallel arrays sorted from best hit to
      # worst hit).  expects an array-like object of hits sorted from best to worst
      # hit with decoys interspersed and a target_setlike object that responds to
      # :include? for the hit object assumes the hit is a decoy if not found
      # in the target set!  if monotonic is false, then the guarantee that
      # qvalues be monotonically increasing is not respected.
      def mixed_target_decoy(best_to_worst, target_setlike, opts={:monotonic => true})
        num_target = 0 ; num_decoy = 0
        monotonic = opts[:monotonic]
        target_hits = []
        qvalues = []
        best_to_worst.each do |hit|
          if target_setlike.include?(hit) 
            num_target += 1
            precision = Ms::ErrorRate::Decoy.precision(num_target, num_decoy)
            target_hits << hit
            qvalues << (1.0 - precision)
          else
            num_decoy += 1
          end
        end
        if opts[:monotonic]
          min_qvalue = qvalues.last 
          qvalues = qvalues.reverse.map do |val| # from worst to best score
            if min_qvalue < val 
              min_qvalue
            else
              min_qvalue = val
              val
            end
          end.reverse
        end
        [target_hits, qvalues]
      end

      module_function :precision

    end
  end
end
