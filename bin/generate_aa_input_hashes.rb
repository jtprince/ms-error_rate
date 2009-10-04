#!/usr/bin/ruby

require 'ms/error_rate/sbv/aa'


if ARGV.size != 3
  puts "usage: #{File.basename(__FILE__)} <min> <aa> <peptide_centric_db>"
  puts "outputs:"
  puts "    <peptide_centric_db>.min<Int>.#{Ms::ErrorRate::Sbv::LENGTH_EXT}"
  puts "    <peptide_centric_db>.min<Int>.#{Ms::ErrorRate::Sbv::AASEQ_EXT}"
  exit
end

(min_num, aa, peptide_centric_db) = ARGV
min_num = min_num.to_i


Ms::ErrorRate::Sbv::Aa.generate_hashes( peptide_centric_db, aa, min_num ).each do |file|
  puts "WROTE: #{file}"
end
