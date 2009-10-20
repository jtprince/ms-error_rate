require 'ms/mascot/dat'

module Ms
  module ErrorRate
    module Qvalue
      module Mascot
      end
    end
  end
end

module Ms::ErrorRate::Qvalue::Mascot
  # returns a hash of arrays of hits (each being: [filename, sequence, charge,
  # mowse, precision]) indexed by the charge state
  # if opts[:z_together], returns a single array
  def qvalues(target_files, decoy_files, opts={})
    opts = {:z_together => false}.merge opts

    all_hits_by_charge = Hash.new {|h,k| h[k] = [] }
    (target_hash, decoy_hash) = [target_files, decoy_files].map do |files|
      hash = {}

      files.each do |file| 
        filename = File.basename(file)

        Ms::Mascot::Dat.open(file) do |dat|
          dat.each_peptide_hit(:yield_nil => false, :with_query => true) do |hit,query|
            hash[hit] = [filename, hit.sequence, query.data['charge'], hit.score] 
            all_hits_by_charge[query.data['charge']] << hit
          end
        end
      end
      hash  # <-- need this for the map
    end

    # Proc.new doesn't do arity checking
    hits_with_qvalues = Proc.new do |hits|
      target_hits_as_arrays_sorted = []
      all_sorted_by_mowse = hits.sort_by{|hit| -(hit.score)}
      (target_hits, qvalues) = Ms::Id::Qvalue.mixed_target_decoy(all_sorted_by_mowse, target_hash)
      target_hits.zip(qvalues) do |hit, qvalue|
        target_hits_as_arrays_sorted.push( target_hash[hit] << qvalue)
      end
      target_hits_as_arrays_sorted
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

  module_function :qvalue


end
