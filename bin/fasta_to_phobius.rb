#!/usr/bin/ruby

require 'mechanize'

page = 'http://phobius.sbc.su.se/'

if ARGV.size == 0
  puts "usage: #{File.basename(__FILE__)} <file>.fasta"
  puts "outputs <file>.phobius "
  puts "in short format"
  puts "submits fasta guys in chunks of #{prot_chunk_size}"
  exit
end


a = WWW::Mechanize.new { |agent|
  agent.user_agent_alias = 'Mac Safari'
}

ARGV.each do |file|
  outfile = file.chomp(File.extname(file)) + '.phobius'
  a.get(page) do |page|
    form = page.forms.first
    form.radiobuttons.select {|v| v.value == 'short' }.first.click
    fu = form.file_uploads.first
    fu.file_name = File.expand_path(file)
    #fu.file_data = IO.read(file)
    reply = form.submit
    html = reply.body
    start = html.index("<pre>") + 5
    stop = html.rindex("</pre>")
    File.open(outfile, 'w') {|out| out.print html[start...stop] }
  end
end

