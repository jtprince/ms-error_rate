#!/usr/bin/ruby -w

require 'rubygems'
require 'optparse'

def require_it(file, gem)
  begin
    require file
  rescue LoadError
    puts "****************************************************"
    puts "  you need to install the '#{gem}' gem:"
    puts "      sudo gem install #{gem}" 
    puts "****************************************************"
    puts $!
    abort "(exiting)"
  end
end

%w(ms/fasta ms/in_silico/digester ms/id).zip(%w(ms-fasta ms-in_silico ms-error_rate)).each do |file, gem|
  require_it file, gem
end

opt = {
  :min_length => 4,
  :missed_cleavages => 2,
  :cleave_initiator_methionine => true,
  :expand_aa => {'X' => Ms::Id::STANDARD_AA},
}
opts = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} <file>.fasta ..."
  op.separator "       returns <file>.msd_clvg<missed_cleavages>.min_aaseq<min_length>.yml"
  op.separator ""
  op.separator "    Initiator Methionines - by default, will generate two peptides"
  op.separator "    for any peptide found at the N-termini starting with 'M'"
  op.separator "    (i.e., one with and one without the leading methionine)"
  op.separator ""
  op.on("--missed-cleavages <Int>", Integer, "max num of missed cleavages (def: #{opt[:missed_cleavages]})") {|v| opt[:missed_cleavages] = v }
  op.on("--min-length <Int>", Integer, "the minimum peptide aaseq length (def: #{opt[:min_length]})") {|v| opt[:min_length] = v }
  op.on("--no-cleaved-methionine", "does not cleave off initiator methionine") { opt[:cleave_initiator_methionine] = false }
  op.on("--no-expand-x", "don't enumerate aa 'X' possibilities") { opt[:expand_aa] = nil }
end

opts.parse!

if ARGV.size == 0
  puts opts.to_s
  exit
end

opt[:remove_digestion_file] = true
opt[:enzyme] = Ms::InSilico::Digester::TRYPSIN

ARGV.each do |file|
  Ms::Id.peptide_centric_db(file, opt)
end
