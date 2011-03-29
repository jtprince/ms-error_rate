#!/usr/bin/env ruby

require 'trollop'
require 'nokogiri'
require 'set'

require 'ms/error_rate/qvalue/pepxml'

opts = Trollop::Parser.new do
  banner %Q{usage: #{File.basename(__FILE__)} <fwd>.xml <decoy>.xml ...
outputs: <fwd>.phq.csv
phq.tsv?: see schema/peptide_hit_qvalues.phq.tsv
}
  opt :z_together, "do not group by charge state", :default => false
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
files.each_slice(2).map do |target, decoy|
  hit_qvalue_pairs = Ms::ErrorRate::Qvalue::Pepxml.target_decoy_qvalues(target, decoy, opt, &:ionscore)
  hits = [] ; qvals = []
  hit_qvalue_pairs.each do |hit, qval|
    hits << hit ; qvals << qval
  end
  outfile = Ms::Ident::PeptideHit::Qvalue.to_phq(target.chomp(File.ext(target)), hits, qvals)
  puts "created: #{outfile}" if $VERBOSE
end


