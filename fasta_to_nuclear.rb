#!/usr/bin/ruby


if ARGV.size == 0
  puts "usage: #{File.basename(__FILE__)} <file>.fasta"
  puts "output: <file>"
  #puts "WARNING!!: you need to run phobius_to_nontransmembrane.rb before"
  #puts "this to weed out transmembrane proteins!"
  exit
end




