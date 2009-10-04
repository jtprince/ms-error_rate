#!/usr/bin/ruby

require 'ms/fasta'
require 'ms/error_rate/sbv/transmembrane'


if ARGV.size != 3
  puts "usage: #{File.basename(__FILE__)} <min_num_certain_tm_seqs> peptide_centric_db phobius_file " 
  puts "outputs:"
  puts "    <peptide_centric_db>.<info>.#{Ms::ErrorRate::Sbv::LENGTH_EXT}"
  puts "    <peptide_centric_db>.<info>.#{Ms::ErrorRate::Sbv::AASEQ_EXT}"
  exit
end

(min_num, peptide_centric_db, phobius) = ARGV
min_num = min_num.to_i

Ms::ErrorRate::Sbv::Transmembrane.generate_hashes( peptide_centric_db, phobius, min_num).each do |file|
  puts "WROTE: #{file}"
end
