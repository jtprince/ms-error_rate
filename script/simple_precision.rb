#!/usr/bin/ruby

require 'optparse'
require 'set'
require 'ms/error_rate/decoy'
require 'ms/mascot/dat'


# returns an array of target hits containing:
# [filename, sequence, charge, mowse, precision]
# sorted by mowse score
def calc_precision(target_files, decoy_files)

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
end

headers = %w(filename aaseq charge mowse precision)

DEF_EXT = "_flip"
NORMAL_EXT = 'precision.yml'
opt = {
  :outfile => NORMAL_EXT,
}

opts = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} <target>.dat <decoy>.dat [<target2>.dat <decoy2>.dat ...]"
  op.separator "for each pair of files"
  op.separator "sorts the peptide hits by mowse score and determines the precision at each hit"
  op.separator ""
  op.separator "writes a yaml file <target>.'#{NORMAL_EXT}' which"
  op.separator "has three keys: 'headers', 'filenames', and 'data'"
  op.separator "    headers contains an array showing what is in the data"
  op.separator "    filenames: (a hash with two keys holding an array of full path names)"
  op.separator "      target:"
  op.separator "      decoy:"
  op.separator "    data: (an array with the data values)"
  op.separator "headers: #{headers.join(', ')}"
  op.separator ""
  op.on("-o", "--outfile <name>", "write to specified file") {|v| opt[:outfile] = v }
  op.on("-a", "--all-together", "process all forwards together and all decoys together", "will output to opt[:outfile] unless -o given") {|v| opt[:all_together] = v }
  op.on("-f", "--find-decoy [ext]", "finds the decoy file, default <file>#{DEF_EXT}.dat", "obviating the need to specify it on the commandline") do |v| 
    if v.is_a? String
      opt[:find_decoy] = v
    else
      opt[:find_decoy] = DEF_EXT
    end
  end
end

opts.parse!

if ARGV.size == 0
  puts opts.to_s
  exit
end


target_files = []
decoy_files = []
if opt[:find_decoy]
  target_files = ARGV.to_a.dup
  decoy_files = target_files.map do |tf|
    ext = File.extname(tf)
    basename = tf.chomp(ext)
    decoy_file = basename + opt[:find_decoy] + ext
    raise ArgumentError, "cannot find #{decoy_file}" unless File.exist?(decoy_file)
    decoy_file
  end
else
  ARGV.each_slice(2) do |target, decoy|
    target_files << target
    decoy_files << decoy
  end
end

filenames = { 'target' => target_files, 'decoy' => decoy_files }

if opt[:all_together]
  target_hits = calc_precision(target_files, decoy_files)
  outfile = opt[:outfile]
  File.open(outfile, 'w') do |out|
    out.print( {'headers' => headers, 'filenames' => filenames, 'data' => target_hits }.to_yaml )
  end
else
  target_files.zip(decoy_files) do |target_file, decoy_file|
    target_hits = calc_precision([target_file], [decoy_file])
    base = target_file.chomp(File.extname(target_file))
    outfile = base + '.' + NORMAL_EXT
    File.open(outfile, 'w') do |out|
      out.print( {'headers' => headers, 'filenames' => filenames, 'data' => target_hits }.to_yaml )
    end
  end
end
