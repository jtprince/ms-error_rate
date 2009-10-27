#!/usr/bin/ruby

require 'orderedhash'
require 'yaml'
require 'set'

if ARGV.size != 2
  puts "usage: #{File.basename(__FILE__)} <gene_ids>.txt summary.yml"
  puts "writes a yml file with unique proteins per qvalue cutoff"
  puts "for each set"
  puts "summary__<setname>__<gene_ids>.yml"
  exit
end

(gene_ids, summary) = ARGV

globs = IO.readlines(gene_ids).reject{|v| v[0,1] == '#'}.map{|v| v.chomp }.select {|v| v =~ /\w/ }

hash = YAML.load_file(summary)
protein_info = hash['protein_info']
results = hash['results']
output_hashes = OrderedHash.new
results.each do |result|

  qvalue_cutoff = result['qvalue_cutoff']
  result['sets'].each do |setname, sethash|
    matches = Set.new
    output_hashes[setname] ||= OrderedHash.new
    proteins = sethash['proteins']
    proteins.each do |ipi,info|
      if info['num_hits_minimal'].first > 0
        all_proteins = [ipi, *info['indistinguishable']]
        all_proteins.each do |id|
          globs.each do |glob|
            if File.fnmatch?(glob, protein_info[id]['Gene_Symbol'])
              matches << protein_info[id]['Gene_Symbol']
            end
          end
        end
      end
    end
    output = matches.to_a.sort
    output_hashes[setname][qvalue_cutoff] = output
  end
end

output_hashes.each do |setname, output|
  gene_ids_base = File.basename(gene_ids, '.*')
  summary_base = summary.chomp(File.extname(summary))
  output_file = [summary_base, setname, gene_ids_base].join("__") + ".yml"

  File.open(output_file, 'w') {|out| out.print output.to_yaml }
end
