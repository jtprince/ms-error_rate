#!/usr/bin/env ruby

if ARGV.size == 0
  puts "usage: #{File.basename(__FILE__)} <file>.fasta ..."
  puts "output: <file>.protid_to_size.yml"
  exit
end

ARGV.each do |file|
  base = file.chomp(File.extname(file))
  outfile = base + ".protid_to_size.yml"
  cnt = 0
  File.open(outfile,'w') do |out|
    Ms::Fasta.foreach(file) do |entry|
      out.puts "#{entry.entry_id}: #{entry.sequence.size}"
      cnt += 1
    end
  end
  puts "wrote #{cnt} protein lengths to: #{outfile}"
end
