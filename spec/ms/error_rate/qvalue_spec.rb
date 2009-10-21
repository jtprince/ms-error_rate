require File.expand_path( File.dirname(__FILE__) + '/../../tap_spec_helper' )

require 'ms/error_rate/qvalue'

Hit = Struct.new(:score, :charge)

class QvalueBaseSpec < MiniTest::Spec

  before(:each) do
    scores = [14,15,13,12,11]
    qvals_expected = [0.5 ,0.0, 2.0/3.0, 3.0/4, 4.0/5]
    @target_hits = scores.zip(Array.new(scores.size, 2)).map {|pair| Hit.new(*pair) } 
    @decoy_hits = scores.zip(Array.new(scores.size, 2)).map {|pair| Hit.new(pair.first-0.5, pair.last) }
    @qval_by_hit = {}
    @target_hits.zip(qvals_expected) {|hit, qval|  @qval_by_hit[hit] = qval }
  end

  it 'can calculate qvalues on target deccoy sets' do
    pairs = Ms::ErrorRate::Qvalue.target_decoy_qvalues(@target_hits, @decoy_hits)
    pairs.each do |hit, qval|
      @qval_by_hit[hit].must_be_close_to(qval, 0.00000001)
    end
  end
  
end
