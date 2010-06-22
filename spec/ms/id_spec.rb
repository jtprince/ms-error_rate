require File.dirname(__FILE__) + '/../spec_helper'

require 'ms/id'

describe 'creating a peptide centric database' do

  it 'can expand out wildcard amino acid combinations (like "X")' do
    array = Ms::Id.expand_peptides('ALXX', 'X' =>  %w(* % &), 'L' => %w(P Q) )
    array.sort.is %w(AP** AP*% AP*& AP%* AP%% AP%& AP&* AP&% AP&& AQ** AQ*% AQ*& AQ%* AQ%% AQ%& AQ&* AQ&% AQ&&).sort
  end
  
end
