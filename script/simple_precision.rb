#!/usr/bin/ruby

require 'optparse'
require 'set'
require 'ms/error_rate/decoy'
require 'ms/mascot/dat'


# returns a hash of arrays of hits (each being: [filename, sequence, charge,
# mowse, precision]) indexed by the charge state
# if opts[:z_together], returns a single array
def mascot_precision(target_files, decoy_files, opts={})
  opts = {:z_together => false}.merge opts

  all_hits_by_charge = Hash.new {|h,k| h[k] = [] }
  (target_hash, decoy_hash) = [target_files, decoy_files].map do |files|
    hash = {}

    files.each do |file| 
      filename = File.basename(file)

      Ms::Mascot::Dat.open(file) do |dat|
        dat.each_peptide_hit(:yield_nil => false, :with_query => true) do |hit,query|
          hash[hit] = [filename, hit.sequence, query.data['charge'], hit.score] 
          all_hits_by_charge[query.data['charge']] << hit
        end
      end
    end
    hash  # <-- need this for the map
  end

  # Proc.new doesn't do arity checking
  hits_with_precision = Proc.new do |hits|
    target_hits_as_arrays_sorted = []
    all_sorted_by_mowse = hits.sort_by{|hit| -(hit.score)}
    (target_hits, precisions) = precision(all_sorted_by_mowse, target_hash)
    target_hits.zip(precisions) do |hit, precision|
      target_hits_as_arrays_sorted.push( target_hash[hit] << precision )
    end
    target_hits_as_arrays_sorted
  end

  if opts[:z_together]
    all_hits = []
    all_hits_by_charge.each do |charge, hits|
      all_hits.push(*hits)
    end
    hits_with_precision.call(all_hits)
  else
    all_hits = []
    all_hits_by_charge.each do |charge,hits|
      all_hits.push(*(hits_with_precision.call(hits)))
    end
    all_hits.sort_by {|v| -(v.last) }
  end
end



# returns [target_hits, precisions] (parallel arrays sorted from best hit to
# worst hit).  expects an array-like object of hits sorted from best to worst
# hit with decoys interspersed and a target_setlike object that responds to
# :include? for the hit object assumes the hit is a decoy if not found in the
# target set!  if monotonic is true, then the precision is guaranteed to be
# the same or increase as score increases, otherwise it is possible to have a
# lower scoring hit with a higher precision since precision is calculated from
# a minimum threshold looking up in score.
def precision(mixed_hits_sorted, target_setlike, opts={:monotonic => true})
  num_target = 0 ; num_decoy = 0
  monotonic = opts[:monotonic]
  target_hits = []
  precisions = []
  mixed_hits_sorted.each do |hit|
    if target_setlike.include?(hit) 
      num_target += 1
      precision = Ms::ErrorRate::Decoy.precision(num_target, num_decoy)
      target_hits << hit
      precisions << precision
    else
      num_decoy += 1
    end
  end
  if opts[:monotonic]
    max_precision = precisions.last 
    precisions = precisions.reverse.map do |val| # from worst to best score
      if max_precision > val 
        max_precision
      else
        max_precision = val
        val
      end
    end.reverse
  end
  [target_hits, precisions]
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
  op.on("--z-together", "combines all charge states for precision calc") {|v| opt[:z_together] = v }
  op.on("-o", "--outfile <name>", "write to specified file") {|v| opt[:outfile] = v }
  op.on("-g", "--group-together", "process all forwards together and all decoys together", "will output to opt[:outfile] unless -o given") {|v| opt[:group_together] = v }
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


print_out = Proc.new do |outfile, filenames, headers, target_hits|
  File.open(outfile, 'w') do |out|
    out.print( {'headers' => headers, 'filenames' => filenames, 'data' => target_hits }.to_yaml )
  end
end

if opt[:group_together]
  filenames = { 'target' => target_files, 'decoy' => decoy_files }
  target_hits= mascot_precision(target_files, decoy_files, opt)
  outfile = opt[:outfile]
  print_out.call(outfile, filenames, headers, target_hits)
else
  target_files.zip(decoy_files) do |target_file, decoy_file|
    filenames = { 'target' => [target_file], 'decoy' => [decoy_file] }
    target_hits = mascot_precision([target_file], [decoy_file], opt)
    base = target_file.chomp(File.extname(target_file))
    outfile = base + '.' + NORMAL_EXT
    print_out.call(outfile, filenames, headers, target_hits)
  end
end
