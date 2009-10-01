#!/usr/bin/ruby

require 'yaml'
require 'set'

extension = ".minimal_prots.yml"

if ARGV.size < 2
  puts "usage: #{File.basename(__FILE__)} <peptide_centric_db> <peptide_seqs> ..."
  puts ""
  puts "outputs the minimal protein set in file '<peptide_seqs>#{extension}"
  puts "<protein>: [<peptide>, <peptide>, <peptide>]"
  puts ""
  puts "peptide_centric_db is in the format: "
  puts "<PEPTIDE>: <ID>-<ID>-<ID>"
  puts ""
  puts "peptide_seqs has one column of peptide sequences"
  puts "(may contain # comments at start of line)"
  exit
end

start = Time.now
pep_db = YAML.load_file(ARGV.shift)
puts "Took #{Time.now - start} seconds to load DB"

ARGV.each do |file|
  basename = file.chomp(File.extname(file))
  outfn = basename + extension
  prots_to_peps = Hash.new {|h,k| h[k] = [] }
  peptides = IO.readlines(file).select {|v| v =~ /^[a-zA-Z]/ }.map {|v| v.chomp}
  peptides_not_found = []
  peptides.each do |seq|
    if pep_db.key?(seq)
      pep_db[seq].split('-').each do |prot|
        prots_to_peps[prot] << seq
      end
    else
      peptides_not_found << seq
    end
  end
  puts "PROTS TO PEPS"
  p prots_to_peps

  puts "not found: "
  p peptides_not_found

  # sort by 1) the total number of proteins and 2) [for ties] the total length
  # of peptides found
  sorted_most_to_least = prots_to_peps.sort_by do |k,v| 
    [ v.size, v.inject(0){|m,s| m+s.size} ]
  end.reverse
  p sorted_most_to_least.size
  p sorted_most_to_least[0,2]

  found_seq = Set.new
  proteins = []

  sorted_most_to_least.each do |prot, peps|
    uniq_peps = peps.select do |pep|
      if found_seq.include?(pep)
        false
      else
        found_seq.add pep
        true
      end
    end
    if uniq_peps.size > 0
      proteins << [prot, uniq_peps]
    end
  end

  File.open(outfn, 'w') do |out|
    proteins.each do |ar|
      out.puts( "#{ar.shift}: [#{ar.join(', ')}]" )
    end
  end
  puts "WROTE: #{outfn}"
end




