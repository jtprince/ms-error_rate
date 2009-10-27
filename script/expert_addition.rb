#!/usr/bin/ruby

require 'yaml'
require 'set'

if ARGV.size == 0
  puts "usage: prog summary__<setname>__name_to_gene_id.yml"
  exit
end

file = ARGV.shift

hash = YAML.load_file(file)

previous_hits = Set.new
results = []
hash.sort.each do |fdr, hits|
  new_hits = hits - previous_hits.to_a
  previous_hits.merge(new_hits)
  results << [fdr, hits.size, *new_hits]
end

results.shift.zip(*results) do |row|
  puts row.join("\t")
end

