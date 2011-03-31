#!/usr/bin/env ruby

require 'trollop'
require 'nokogiri'
require 'set'

require 'ms/ident/pepxml'
require 'ms/ident/peptide_hit/qvalue'

EXT = Ms::Ident::PeptideHit::Qvalue::FILE_EXTENSION
combine_base  = "combined"

opts = Trollop::Parser.new do
  banner %Q{usage: #{File.basename(__FILE__)} <target>.xml <decoy>.xml ...
outputs: <fwd>.phq.tsv
phq.tsv?: see schema/peptide_hit_qvalues.phq.tsv
}
  opt :combine, "groups target and decoy hits together from all files, writing to #{combine_base}#{EXT}", :default => false
  opt :z_together, "do not group by charge state", :default => false
  opt :decoy_first, "decoy files are before target files", :default => false
  opt :verbose, "be verbose", :default => false
end

opt = opts.parse(ARGV)
if ARGV.size == 0 || (ARGV.size%2 != 0)
  puts "\n\n!! only even numbers of files accepted (target decoy target decoy ...) !!\n\n" if (ARGV.size%2 != 0)
  opts.educate
  exit 
end

$VERBOSE = opt.delete(:verbose)

files = ARGV.to_a

files.each_slice(2).map!(&:reverse) if opt[:decoy_first]

groups_of_search_hits = files.map do |file|
  Ms::Ident::Pepxml.search_hits(file)
end

to_run = {}
if opt[:combine]
  all_target = [] ; all_decoy = []
  groups_of_search_hits.each_slice(2) do |target_hits,decoy_hits| 
    all_target.push(*target_hits) ; all_decoy.push(*decoy_hits)
  end
  to_run[combine_base] = [all_target, all_decoy]
else
  groups_of_search_hits.each_slice(2).zip(files) do |target_hits_and_decoy_hits, file|
    to_run[file.chomp(File.extname(file)) + EXT] = target_hits_and_decoy_hits
  end
end

to_run.each do |file, target_decoy_pair|
  hit_qvalue_pairs = Ms::ErrorRate::Qvalue.target_decoy_qvalues(target_decoy_pair.first, target_decoy_pair.last, :z_together => opt[:z_together])
  hits = [] ; qvals = []
  hit_qvalue_pairs.each do |hit, qval|
    hits << hit ; qvals << qval
  end
  outfile = Ms::Ident::PeptideHit::Qvalue.to_file(file, hits, qvals)
  puts "created: #{outfile}" if $VERBOSE
end


