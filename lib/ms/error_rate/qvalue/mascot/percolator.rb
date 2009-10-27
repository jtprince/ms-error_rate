
require 'ms/mascot/dat'
require 'ms/error_rate/qvalue'
require 'ms/error_rate/qvalue/mascot'

module Ms
  module ErrorRate
    module Qvalue
      module Mascot
        module Percolator

          module_function
          # returns an array of Structs where the keys are the first line
          # everything is cast properly
          # three additional keys are available query_num, rank, sequence
          # sequence is the amino acid sequence without the surrounding X's
          # and dots.
          # (with '-' substituted for '_')
          def tab_txt(file)
            hits = []
            File.open(file) do |io|
              # PSMId	score	q-value	posterior_error_prob	peptide	proteinIds
              atts = io.gets.chomp.split("\t").map {|v| v.gsub('-', '_').to_sym } 
              atts.push(:query_num, :rank, :sequence)
              struct_class = Struct.new("Hit", *atts)
              
              io.each do |line|
                (query_rank, score, qvalue, perrp, peptide, *prots ) = line.chomp.split("\t")
                (query, rank) = query_rank.split(';').map {|v| v.split(':').last.to_i }
                
                hits << struct_class.new(query_rank, score.to_f, qvalue.to_f, perrp.to_f, peptide, prots, query, rank, peptide.split('.')[1])
              end
            end
            hits
          end

        end
      end
    end
  end
end

module Ms::ErrorRate::Qvalue::Mascot::Percolator

  module_function
  # returns an array of Structs of PeptideHit(:filename, :query_title, :charge, :sequence, :mowse, :qvalue)
  # opts = 
  #   :min_peptide_length => Integer
  def qvalues(datp_files, tab_txt_files, opts={})
    min_pep_len = opts[:min_peptide_length]

    # we only want the top hit per query title (which should ensure that we
    # get the top hit per scan)
    hits_by_query_title = Hash.new {|h,k| h[k] = [] }
    datp_files.zip(tab_txt_files) do |datp_file, tab_txt_file| 
      # build a hash based on the sequence
      structs = Ms::ErrorRate::Qvalue::Mascot::Percolator.tab_txt( tab_txt_file )
      qvalue_by_query_rank = {}
      structs.each do |struct|
        qvalue_by_query_rank[[struct.query_num, struct.rank]] = struct.q_value
      end

      base_no_ext = File.basename(datp_file, '.*')
      Ms::Mascot::Dat.open(datp_file) do |dat|
        dat.each_peptide_hit(:by => :groups, :yield_nil => false, :with_query => true) do |hits,query|
          hits.each do |hit|
            if qval = qvalue_by_query_rank[[hit.query_num, hit.hit_num]]
              hit_as_struct = Ms::ErrorRate::Qvalue::Mascot::MascotPeptideHit.new(base_no_ext, query.title, query.charge, hit.sequence, hit.score, qval)
              hits_by_query_title[hit_as_struct.query_title] << hit_as_struct
            end
          end
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
end
