#!/usr/bin/ruby

require 'optparse'

require 'ms/error_rate/sbv/peptide_based'
require 'ms/error_rate/sbv/protein_based'

opts = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} peptide_centric_db [OPTION]" 
  op.on("--tm <phobius,min>", Array, "transmembrane, <phobius> is path to phobius ", "output file (see fasta_to_phobius.rb)" "<min> is the min number of tm sequences required") {|v| opt[:tm] = [v.first, v.last.to_i]}
  op.on("--aa <aa,min>", Array, "amino acid, <aa> is a string found in the peptides", "<min> is the min number of required for counting") {|v| opt[:aa] = [v.first, v.last.to_i]}
  op.on("--protein-bias <name,file>", Array, "<name> bias, <file> is path to a yaml hash", "    keyed prot -> <0-1>") {|v| opt[v.first.to_sym] = [v.last]}
  puts "outputs for each bias type:"
  puts "    <peptide_centric_db>.<info>.#{Ms::ErrorRate::Sbv::LENGTH_EXT}"
  puts "    <peptide_centric_db>.<info>.#{Ms::ErrorRate::Sbv::AASEQ_EXT}"
end

opts.parse!

if ARGV.size == 0
  puts opts.to_s
  exit
end


klass = Ms::ErrorRate::Sbv
tm_klass = Ms::ErrorRate::Sbv::Transmembrane

if op[:tm]
  index = TransmembraneIndex.new(op[:tm].first)

  protid_to_transmembrane = {}
  regexp = nil
  index.each do |k,v|
    regexp ||= Ms::Fasta.id_regexp(k)
    new_key = regexp.match(k)[1] 
    protid_to_transmembrane[new_key] = ((v[:num_certain_transmembrane_segments] >= op[:tm].last) ? 1 : 0)
  end

  fnames = tm_klass.generate_hashes( peptide_centric_db, protid_to_transmembrane)
  end
  fnames.each do |file|
    puts "WROTE: #{file}"
  end
end


if op[:aa]
  klass.generate_hashes( peptide_centric_db, 

end



