#!/usr/bin/ruby

require 'ms/fasta'
require 'transmembrane/phobius.rb'

if ARGV.size != 3
  puts "usage: #{File.basename(__FILE__)} max_num_tm phobius_short_file <file>.fasta"
  puts "max_num_tm = max # of transmembrane sequences allowed to be a non-transmembrane."
  puts ""
  puts "outputs: <file>_NONTM.fasta"
end

(max_num_tm, phobius_short_file, fasta_db_file) = ARGV
max_num_tm = max_num_tm.to_i

base = fasta_db_file.chomp(File.extname(fasta_db_file))
outfile = base + "_NONTM.fasta"

index = Phobius::Index.new(phobius_short_file)

File.open(outfile, 'w') do |out|
  Ms::Fasta.open(fasta_db_file) do |fasta|
    fasta.each do |entry|
      key = index.reference_to_key(entry.header)
      abort "can't find key: #{key} for #{entry.header}"  unless index.key?(key)
      num_tms = index[key][:num_certain_transmembrane_segments]
      if num_tms <= max_num_tm
        out.print entry.to_s
      end
    end
  end
end


