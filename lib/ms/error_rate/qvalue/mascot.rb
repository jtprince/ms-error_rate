require 'ms/error_rate/qvalue'
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
  MEMBERS = [:filename, :query_title, :charge, :sequence, :mowse, :qvalue]
  MascotPeptideHit = Struct.new(*MEMBERS) do 
    # emits an array rather than a Struct object
    def to_yaml(*args)
      to_a.to_yaml(*args)
    end
  end

  module_function
  # returns an array of Structs of PeptideHit(:filename, :query_title, :charge, :sequence, :mowse, :qvalue)
  # opts = 
  #   :min_peptide_length => Integer
  def qvalues(target_files, decoy_files, opts={})
    min_pep_len = opts[:min_peptide_length]

    # we only want the top hit per query title (which should ensure that we
    # get the top hit per scan)
    (target_hits, decoy_hits) = [target_files, decoy_files].map do |files|
      hits_by_query_title = Hash.new {|h,k| h[k] = [] }
      files.each do |file| 
        base_no_ext = File.basename(file, '.*')
        Ms::Mascot::Dat.open(file) do |dat|
          dat.each_peptide_hit(:by => :top, :yield_nil => false, :with_query => true) do |hit,query|

            hit_as_struct = MascotPeptideHit.new(base_no_ext, query.title, query.charge, hit.sequence, hit.score)
            hits_by_query_title[hit_as_struct.query_title] << hit_as_struct
          end
        end
      end

      final_hits = []
      hits_by_query_title.each do |title, hits|
        best_hit = 
          if hits.size == 1
            hits.first
          else
            hits.sort_by(&:mowse).last
          end
        # FILTER HERE:
        # ONLY TAKE the BEST HIT IF it passes any filters
        if min_pep_len
          next unless best_hit.sequence.size >= min_pep_len
        end
        final_hits << best_hit
      end
      final_hits
    end
    pairs = Ms::ErrorRate::Qvalue.target_decoy_qvalues(target_hits, decoy_hits, opts, &:mowse)
    pairs.map do |hit, qval| 
      hit.qvalue = qval 
      hit
    end
  end
end
