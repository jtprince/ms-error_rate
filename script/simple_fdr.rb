#!/usr/bin/ruby

require 'set'
require 'ms/error_rate/decoy'
require 'ms/mascot/dat'

if ARGV.size < 3
  puts "usage: #{File.basename(__FILE__)} <fdr_cutoff> <target>.dat <decoy>.dat [<target2>.dat <decoy2>.dat ...]"
  puts "sorts the peptide hits by mowse score and determines the number of unique"
  puts "peptide sequences and unique peptides at a given charge state before cutoff"
  puts "writes these to a file"
  exit
end

cutoff = ARGV.shift.to_f

target_files = []
decoy_files = []
ARGV.each_slice(2) do |target, decoy|
  target_files << target
  decoy_files << decoy
end

if ARGV.size > 3
  puts "TARGET: "
  target_files.each {|v| puts v }
  puts "DECOY: "
  decoy_files.each {|v| puts v }
end

all_hits = []
(target_hash, decoy_hash) = [target_files, decoy_files].map do |files|
  hash = {}

  files.each do |file| 
    Ms::Mascot::Dat.open(file) do |dat|
      dat.each_peptide_hit(:yield_nil => false, :with_query => true) do |hit,query|
        hash[hit] = query.data['charge']
        all_hits << hit
      end
    end
  end
  hash
end

sorted_by_highest_mowse = all_hits.sort_by {|v| v.score }.reverse

num_target = 0
num_decoy = 0

most_hits = sorted_by_highest_mowse.map do |hit|
  if target_hash.key?(hit) 
    num_target += 1
    precision = Ms::ErrorRate::Decoy.precision(num_target, num_decoy)
    p precision
    [precision, hit, target_hash[hit]]
  elsif decoy_hash.key?(hit)
    num_decoy += 1
    nil
  else
    raise "hit must be target or decoy I reckon"
  end
end.compact.reverse

(below, above) =  most_hits.partition do |ar|
  ar.first < cutoff
end

unique_sequences = Set.new
unique_ions = Set.new
above.each do |ar|
  sequence = ar[1].sequence
  seq_plus_charge = sequence + ar.last
  unique_sequences.add sequence
  unique_ions.add  seq_plus_charge
end

puts "at a #{(1.0-cutoff)*100} % FDR: "
puts "  num psms: #{above.size}"
puts "  num unique seqs: #{unique_sequences.size}"
puts "  num unique ions: #{unique_ions.size}"



