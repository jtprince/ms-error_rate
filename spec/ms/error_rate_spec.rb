require File.dirname(__FILE__) + '/../spec_helper'

require 'ms/error_rate'
require 'ostruct'

xdescribe 'not quite sure what this is' do

  it 'calculates bayesian probabilities' do
    # C = is a correct ID
    # T = transmembrane content
    # Y = cysteine content
    # A = abundance
    # p(C|T,Y,A) = p(T|C)p(Y|C)p(A|C)p(C) / p(T)p(Y)p(A)
    peps.map do |pep|
      # what is the probability of that un-transmembraneyness being correct?
      # what is the probability of that un-cysteineness being correct?
      # what is the probability of that high abundanceness being correct?
      pep.bayes_probs.reduce(prob_being_correct) do |prob| 
      end
      p_correct = pep.prior_prob_correct
      pep.not_transmembrane? * pep.not_cysteine? * pep.not_low_abundance?
    end
  end
  
end
