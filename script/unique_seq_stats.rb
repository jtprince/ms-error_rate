#!/usr/bin/ruby -w

require 'yaml'
require 'optparse'

opt = {}
opts = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} <precision_file>.yml ..."
  op.separator "outputs yaml formatted information:"
  op.separator "---"
  op.separator "filename: <pathgiven>"
  op.separator "num unique aaseqs: <Int>"
  op.separator "num unique aaseqs+charge: <Int>"
  op.separator "num peptide hits: <Int>"
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

ARGV.each do |file|
  hash = YAML.load_file(file)

  prec_index = hash['headers'].index('precision')
  mowse_index = hash['headers'].index('mowse')

  hits = hash['data']
  low_to_hi_by_mowse = hits.sort_by {|v| v[mowse_index] }

  prec_cutoff_index = nil
  low_to_hi_by_mowse.each_with_index do |hit,i|
    if hit[prec_index] >= opt[:cutoff]
      prec_cutoff_index = i
      break
    end
  end

  above_cutoff = low_to_hi_by_mowse[prec_cutoff_index..-1].reverse

  unique_sequences = Set.new
  unique_ions = Set.new
  above.each do |ar|
    sequence = ar[1].sequence
    seq_plus_charge = sequence + ar.last
    unique_sequences.add sequence
    unique_ions.add  seq_plus_charge
  end

  puts "at a #{precision_cutoff} precision (#{(1.0-precision_cutoff)*100} % FDR): "
  puts "  num psms: #{above.size}"
  puts "  num unique seqs: #{unique_sequences.size}"
  puts "  num unique ions: #{unique_ions.size}"



  File.open(seqfile_full, 'w') do |out|
    out.print "# target files: "
    out.puts target_files.join(',')

  out.print "# decoy files: "
  out.puts decoy_files.join(',')

  unique_sequences.to_a.sort.each do |seq|
    out.puts seq
  end
end
puts "wrote sequences to: #{seqfile_full}"
