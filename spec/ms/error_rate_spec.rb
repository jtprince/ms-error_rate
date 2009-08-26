require File.expand_path( File.dirname(__FILE__) + '/../tap_spec_helper' )

require 'ms/error_rate'
require 'ostruct'

class ErrorRateBase < MiniTest::Spec

  before(:each) do
  end

  it 'calculates bayesian probabilities' do
    # p(a|b) = p(b|a)p(a) / p(b)
    p_b_given_a = 0.6 # prob of sample info given its correct
    peps.map do |pep|
      p_correct = pep.prior_prob_correct
      pep.not_transmembrane? * pep.not_cysteine? * pep.not_low_abundance?
    end
  end
  
end
