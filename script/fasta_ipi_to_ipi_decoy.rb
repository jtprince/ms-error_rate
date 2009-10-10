#!/usr/bin/ruby

if ARGV.size == 0
  puts "usage: #{File.basename(__FILE__)} <IPI_based>.fasta ..."
  puts "moves any leading \"><.*_>\" to the IPI value"
  puts "for example:"
  puts ">DCY_IPI:IPI0032311.1|STUFF ->  >IPI:DCY_IPI0032311.1|STUFF"
  exit
end

ARGV.each do |file|
  tmp = file + '.tmp'
  if File.exist?(tmp) ; warn "Skipping #{file} since #{tmp} exists" ; next  end
  File.open(tmp, 'w') do |out|
    IO.foreach(file) do |line|
      if line =~ />([^\:\|]+_)/
        line.sub!("#{$1}IPI:IPI", "IPI:#{$1}IPI")
      end
      out.print line
    end
  end
  FileUtils.mv tmp, file
end
