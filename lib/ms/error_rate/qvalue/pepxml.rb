require 'ms/error_rate/qvalue'

module Ms ; end
module Ms::ErrorRate ; end
module Ms::ErrorRate::Qvalue ; end

module Ms::ErrorRate::Qvalue::Pepxml
  DELIMITER = "\t"
  PeptideHit = Struct.new(:aaseq, :charge, :ionscore, :qvalue)

  module_function

  # writes a phq.tsv file based on the target's basename
  # retrieves the aaseq, charge, and all search_score keys and values for use
  # in the search_hit.  caller must provide a sort_by block, where the best
  # hits are last.  charge is an integer, and all other search scores are cast
  # as floats.  returns the output filename.
  def to_phq(target_pepxml, decoy_pepxml, opt={}, &sort_by)

    # this is a list of high quality peptide hits associated with each group
    fields = [:aaseq, :charge]
    ss_names = []
    have_ss_names = false
    (target_hits, decoy_hits) = [target_pepxml, decoy_pepxml].map do |file|
      # begin with aaseq, charge
      File.open(file) do |io|
        doc = Nokogiri::XML.parse(io, nil, nil, Nokogiri::XML::ParseOptions::DEFAULT_XML | Nokogiri::XML::ParseOptions::NOBLANKS)
        # we can work with namespaces, or just remove them ...
        doc.remove_namespaces!
        root = doc.root
        search_hits = root.xpath('//search_hit')
        search_hits.map do |search_hit| 
          aaseq = search_hit['peptide']
          charge = search_hit.parent.parent['assumed_charge'].to_i
          search_score_nodes = search_hit.children.select {|node| node.name == 'search_score' }
          ss_values = []
          search_score_nodes.each do |node|
            ss_names << node['name'].to_sym unless have_ss_names
            ss_values << node['value'].to_f
          end
          have_ss_names = true
          [aaseq, charge, *ss_values]
        end
      end
    end

    fields.push(*ss_names)
    peptide_hit_class = Struct.new(*fields)
    (t_hits, d_hits) = [target_hits, decoy_hits].map {|hits| hits.map {|hit_values| peptide_hit_class.new(*hit_values) } }

    hit_qvalue_pairs = Ms::ErrorRate::Qvalue.target_decoy_qvalues(t_hits, d_hits, :z_together => opt[:z_together], &sort_by)

    newfile = target_pepxml.chomp(File.extname(target_pepxml)) + ".phq.tsv"
    File.open(newfile,'w') do |out| 
      out.puts %w(aaseq charge qvalue).join(DELIMITER)
      hit_qvalue_pairs.each do |hit, qvalue|
        out.puts [hit.aaseq, hit.charge, qvalue].join(DELIMITER)
      end
    end
    newfile
  end
end
