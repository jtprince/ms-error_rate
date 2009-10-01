#!/usr/bin/ruby -w

require 'set'
require 'yaml'
require 'optparse'

opt = {}
opts = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} <precision_file>.yml ..."
  op.separator "outputs information collected by combining hits from files:"
  op.separator "---"
  op.separator "filenames: "
  op.separator "- <pathgiven>"
  op.separator "num_unique_aaseqs: <Int>"
  op.separator "num_unique_aaseqs_charge: <Int>"
  op.separator "num_peptide_hits: <Int>"
  op.separator ""
  op.separator "NOTE: if a precision cutoff is given, all hits that have a better"
  op.separator "score than the worst score at the cutoff are included, even if "
  op.separator "the precision for that hit was below the cutoff"
  op.separator "this prevents early, local aberrations in precision from messing"
  op.separator "up the analysis"
  op.separator ""
  op.on("-p", "--precision <0-1>", Float, "precision cutoff") {|v| opt[:cutoff] = v }
  op.on("-f", "--fdr <0-1>", Float, "false discovery rate cutoff (1-precision)") {|v| opt[:cutoff] = 1.0 - v }
end

opts.parse!

if ARGV.size == 0
  puts opts.to_s  
  exit
end

unique_sequences = Set.new
unique_ions = Set.new
all_hits = []

ARGV.each do |file|
  hash = YAML.load_file(file)

  prec_index = hash['headers'].index('precision')
  mowse_index = hash['headers'].index('mowse')
  aaseq_index = hash['headers'].index('aaseq')
  charge_index = hash['headers'].index('charge')

  above_cutoff.each do |ar|
    sequence = ar[aaseq_index]
    seq_plus_charge = sequence + ar[charge_index]
    unique_sequences.add sequence
    unique_ions.add  seq_plus_charge
  end
end

prec_k = 'precision cutoff'
fn_k = 'filenames'
uniq_aaseq_k = 'num unique aaseqs'
uniq_ions_k = 'num unique aaseqs+charge'
num_hits_k = 'num peptide hits'

order = [fn_k, prec_k, num_hits_k, uniq_ions_k, uniq_aaseq_k]

results = {}
results[fn_k] = '[' + ARGV.join(", ") + ']'
results[prec_k] = opt[:cutoff]
results[uniq_aaseq_k] = unique_sequences.size
results[uniq_ions_k] = unique_ions.size
results[num_hits_k] = all_hits.size

order.each do |key|
  puts "#{key}: #{results[key]}"
end
