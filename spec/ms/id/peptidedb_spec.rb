require File.dirname(__FILE__) + '/../../spec_helper'

require 'yaml'
path = 'ms/id/peptidedb'
require path 

module Kernel
 
  def capture_stdout
    out = StringIO.new
    $stdout = out
    yield
    out.rewind
    return out.read
  ensure
    $stdout = STDOUT
  end
 
end

describe 'amino acid expansion' do

  it 'can expand out wildcard amino acid combinations' do
    array = Ms::Id::Peptidedb.expand_peptides('ALXX', 'X' =>  %w(* % &), 'L' => %w(P Q) )
    array.sort.is %w(AP** AP*% AP*& AP%* AP%% AP%& AP&* AP&% AP&& AQ** AQ*% AQ*& AQ%* AQ%% AQ%& AQ&* AQ&% AQ&&).sort
  end

  it 'will not expand explosive combinations (>MAX_NUM_AA_EXPANSION)' do
    # this is from real data
    worst_case = 'LTLLRPEKHEAATGVDTICTHRVDPIGPGLXXEXLYWELSXLTXXIXELGPYTLDR'
    Ms::Id::Peptidedb.expand_peptides(worst_case, 'X' =>  %w(* % &)).nil?.is true
  end

  it 'returns the peptide in the array if no expansion' do
    array = Ms::Id::Peptidedb.expand_peptides('ZZZZZ', 'X' =>  %w(* % &), 'L' => %w(P Q) )
    array.is ['ZZZZZ']
  end

end

describe 'creating a peptide centric database' do

  before do
    @fasta_file = [TESTFILES, path, 'uni_11_sp_tr.fasta'].join('/')
    #@output_file = [TESTFILES, path, 'uni_11_sp_tr.'].join('/')
    @output_file = [TESTFILES, path, "uni_11_sp_tr.msd_clvg2.min_aaseq4.yml"].join('/')
  end

  it 'converts a fasta file into peptide centric db' do
    Ms::Id::Peptidedb.cmdline([@fasta_file])
    ok File.exist?(@output_file)
    sorted = YAML.load_file(@output_file).sort
    # these are merely frozen, not perfectly defined
    sorted.first.is ["AAFDDAIAELDTLSEESYK", ["P62258"]]
    sorted.last.is ["YWCRLGPPRWICQTIVSTNQYTHHR", ["D2KTA8"]]
    sorted.size.is 728
    File.unlink(@output_file)
  end

  it 'lists approved enzymes and exits' do
    output = capture_stdout do
      begin
        Ms::Id::Peptidedb.cmdline(['--list-enzymes'])
      rescue SystemExit
        1.is 1 # we exited
      end
    end
    lines = output.split("\n")
    ok lines.include?("trypsin")
    ok lines.include?("chymotrypsin")
  end
end
