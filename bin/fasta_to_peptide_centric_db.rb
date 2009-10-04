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
  :min_length => 8,
  :missed_cleavages => 1,
}
opts = OptionParser.new do |op|
  op.banner = "usage: #{File.basename(__FILE__)} missed_cleavages <file>.fasta ..."
  op.separator "       returns <file>.msd_clvg<missed_cleavages>.min_aaseq<min_length>.yml"
  op.separator ""
  op.on("--missed-cleavages <Int>", Integer, "max num of missed cleavages (def: #{opt[:missed_cleavages]})") {|v| opt[:missed_cleavages] = v }
  op.on("--min-length <Int>", Integer, "the minimum peptide aaseq length (def: #{opt[:min_length]})") {|v| opt[:min_length] = v }
  
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
