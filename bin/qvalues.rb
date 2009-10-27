#!/usr/bin/ruby

require 'optparse'
require 'ms/error_rate/qvalue'

DEF_EXT = "_flip"
NORMAL_EXT = 'qval.yml'

def print_out(outfile, filenames, headers, target_hits)
  File.open(outfile, 'w') do |out|
    out.print( {'headers' => headers, 'filenames' => filenames, 'data' => target_hits }.to_yaml )
  end
end

opt = {
  :outfile => NORMAL_EXT,
  :min_peptide_length => 9,
}

opts = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} <target> <decoy> [... (as pairs)]"
  op.separator "or:    #{File.basename(__FILE__)} <target>.datp <target>.tab.txt [... (as pairs)]"
  op.separator "for each pair of files"
  op.separator "sorts the peptide hits by score and determines the precision at each hit"
  op.separator ""
  op.separator "writes a yaml file <target>.'#{NORMAL_EXT}' which"
  op.separator "has three keys: 'headers', 'filenames', and 'data'"
  op.separator "    headers contains an array showing what is in the data"
  op.separator "    filenames: (a hash with two keys holding an array of full path names)"
  op.separator "      target:"
  op.separator "      decoy:"
  op.separator "    data: (an array with the data values)"
  op.separator "headers:  <the headers of the hits>"
  op.separator ""
  op.separator "headers guaranteed to have at least: filename, query_title, charge, sequence, qvalue"
  op.separator ""
  op.on("-l", "--min-peptide-length <Int>", Integer, "min num aa's to accept (default: #{opt[:min_peptide_length]})") {|v| opt[:min_peptide_length] = v }
  op.on("--z-together", "combines all charge states for precision calc") {|v| opt[:z_together] = v }
  op.on("-o", "--outfile <name>", "write to specified file") {|v| opt[:outfile] = v }
  op.on("-g", "--group-together", "process all forwards together and all decoys together", "will output to opt[:outfile] unless -o given") {|v| opt[:group_together] = v }
  op.on("-f", "--find-decoy [ext]", "finds the decoy file, default <file>#{DEF_EXT}.<ext>", "obviating the need to specify it on the commandline") do |v| 
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

require 'ms/error_rate/qvalue/mascot'
require 'ms/error_rate/qvalue/mascot/percolator'

mascot_percolator = (File.extname(target_files.first) == '.datp')
headers = Ms::ErrorRate::Qvalue::Mascot::MEMBERS.map(&:to_s)
if opt[:group_together]
  outfile = opt[:outfile]
  if mascot_percolator
    filenames = { 'target' => target_files, 'decoy' => decoy_files }
    # in the case of mascot_percolator, the "target" files are .datp files and
    # "decoy" files the .tab.txt files
    target_hits = Ms::ErrorRate::Qvalue::Mascot::Percolator.qvalues( target_files, decoy_files, opt).sort_by(&:qvalue)
  else
    filenames = { 'target' => target_files, 'decoy' => decoy_files }
    target_hits = Ms::ErrorRate::Qvalue::Mascot.qvalues(target_files, decoy_files, opt).sort_by(&:qvalue)
  end
  print_out(outfile, filenames, headers, target_hits)
else
  target_files.zip(decoy_files) do |target_file, decoy_file|
    if mascot_percolator
      filenames = { 'datp' => [target_file], 'tab_txt' => [decoy_file] }
      target_hits = Ms::ErrorRate::Qvalue::Mascot::Percolator.qvalues([target_file], [decoy_file], opt).sort_by(&:qvalue)
    else
      filenames = { 'target' => [target_file], 'decoy' => [decoy_file] }
      target_hits = Ms::ErrorRate::Qvalue::Mascot.qvalues([target_file], [decoy_file], opt).sort_by(&:qvalue)
    end
    base = target_file.chomp(File.extname(target_file))
    outfile = base + '.' + NORMAL_EXT
    print_out(outfile, filenames, headers, target_hits)
  end
end
