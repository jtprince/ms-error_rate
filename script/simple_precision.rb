#!/usr/bin/ruby

require 'optparse'
require 'set'
require 'ms/error_rate/decoy'
require 'ms/mascot/dat'

headers = %w(filename aaseq charge mowse precision)

NORMAL_EXT = 'precision.yml'
opt = {
  :outfile => NORMAL_EXT,
}

opts = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} <target>.dat <decoy>.dat [<target2>.dat <decoy2>.dat ...]"
  op.separator "sorts the peptide hits by mowse score and determines the precision at each hit"
  op.separator ""
  op.separator "writes a yaml file '#{opt[:outfile]}'"
  op.separator "has three keys, 'headers', 'filenames', and 'data'"
  op.separator "    headers contains an array showing what is in the data"
  op.separator "    filenames: (a hash with two keys holding an array of full path names)"
  op.separator "      target:"
  op.separator "      decoy:"
  op.separator "    data: (an array with the data values)"
  op.separator "headers: #{headers.join(', ')}"
  op.separator ""
  op.on("-o", "--outfile <name>", "write to specified file") {|v| opt[:outfile] = v }
  op.on("-m", "--mimic", "outfile is basename of first file + '.#{NORMAL_EXT}'") {|v| opt[:mimic] = v }
end

opts.parse!

if ARGV.size < 2
  puts opts.to_s
  exit
end


target_files = []
decoy_files = []
ARGV.each_slice(2) do |target, decoy|
  target_files << target
  decoy_files << decoy
end

filenames = { 'target' => target_files, 'decoy' => decoy_files }

all_hits = []
(target_hash, decoy_hash) = [target_files, decoy_files].map do |files|
  hash = {}

  files.each do |file| 
    filename = File.basename(file)

    Ms::Mascot::Dat.open(file) do |dat|
      dat.each_peptide_hit(:yield_nil => false, :with_query => true) do |hit,query|
        hash[hit] = [filename, hit.sequence, query.data['charge'], hit.score] 
        all_hits << hit
      end
    end
  end
  hash
end

sorted_by_highest_mowse = all_hits.sort_by {|v| v.score }.reverse

num_target = 0
num_decoy = 0

target_hits = sorted_by_highest_mowse.map do |hit|
  if target_hash.key?(hit) 
    num_target += 1
    precision = Ms::ErrorRate::Decoy.precision(num_target, num_decoy)
    target_hash[hit] << precision
  elsif decoy_hash.key?(hit)
    num_decoy += 1
    nil
  else
    raise "hit must be target or decoy I reckon"
  end
end.compact

outfile = 
  if opt[:mimic]
    ftf = target_files.first
    base = ftf.chomp(File.extname(ftf))
    base + '.' + NORMAL_EXT
  else
    opt[:outfile]
  end

File.open(outfile, 'w') do |out|
  out.print( {'headers' => headers, 'filenames' => filenames, 'data' => target_hits }.to_yaml )
end
