#!/usr/bin/env ruby

require 'trollop'
require 'nokogiri'
require 'set'

require 'ms/error_rate/qvalue'

opts = Trollop::Parser.new do
  banner %Q{usage: #{File.basename(__FILE__)} <fwd>.xml <decoy>.xml ...
outputs: <fwd>.phq.csv
phq.tsv?: see schema/peptide_hit_qvalues.phq.tsv
}
  opt :z_together, "do not group by charge state", :default => false
end

DELIMITER = "\t"

opt = opts.parse(ARGV)
if ARGV.size == 0 || (ARGV.size%2 != 0)
  puts "\n\n!! only even numbers of files accepted (target decoy target decoy ...) !!\n\n" if (ARGV.size%2 != 0)
  opts.educate
  exit 
end

files = ARGV.to_a

PeptideHit = Struct.new(:aaseq, :charge, :ionscore, :qvalue)

# this is a list of high quality peptide hits associated with each group
peptide_hits_per_file = files.map do |file|
  File.open(file) do |io|
    doc = Nokogiri::XML.parse(io, nil, nil, Nokogiri::XML::ParseOptions::DEFAULT_XML | Nokogiri::XML::ParseOptions::NOBLANKS)
    # we can work with namespaces, or just remove them ...
    doc.remove_namespaces!
    root = doc.root
    search_hits = root.xpath('//search_hit')
    search_hits.map do |search_hit| 
      aaseq = search_hit['peptide']
      ionscore = search_hit.children.find {|node| node.name == 'search_score' && node['name'] == 'ionscore' }['value'].to_f
      charge = search_hit.parent.parent['assumed_charge'].to_i
      PeptideHit.new(aaseq, charge, ionscore)
    end
  end
end

hits_per_target = peptide_hits_per_file.each_slice(2).map do |target_hits, decoy_hits|
  pairs = Ms::ErrorRate::Qvalue.target_decoy_qvalues(target_hits, decoy_hits, :z_together => opt[:z_together], &:ionscore)
  target_peptide_hits = pairs.map {|peptide_hit, qvalue| peptide_hit.qvalue = qvalue ; peptide_hit }
end

files.each_slice(2).map(&:first).zip(hits_per_target) do |file, hits|
  newfile = file.chomp(File.extname(file)) + ".phq.tsv"
  File.open(newfile,'w') do |out| 
    out.puts %w(aaseq charge qvalue).join(DELIMITER)
    hits.each do |hit|
      out.puts hit.values_at(0,1,3).join(DELIMITER)
    end
  end
end

