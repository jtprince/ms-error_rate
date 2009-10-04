require 'set'

class Reverser
  attr_accessor :obj

  def initialize(obj)
    @obj = obj
  end

  def <=>(other)
    other.obj <=> self.obj
  end
end

class Object
  def rev
    Reverser.new(self)
  end
end

module Enumerable
  # Provides sorting on multiple attributes (each directional) where atts is 
  # an array of symbols.
  # the default is to sort ascending (small to large).
  # the option :down => Symbol or ArrayOfSymbols
  #   sort_by_attributes(:age,:height,:weight) # -> sorts by age, height, and weight
  #   sort_by_attributes(:age,:height,:weight, :down => :height) # -> same as above, but sorts height from large to small
  #   sort_by_attributes(:age,:height,:weight, :down => [:height,:weight]) # -> same as above, but sorts height and weight from large to small
  def sort_by_attributes(*atts)
    down = 
      if atts.last.is_a? Hash
        hash = atts.pop
        unless hash[:down].is_a?(Array)
          hash[:down] = [hash[:down]]
        end
        Set.new(hash[:down])
      else
        Set.new
      end
    self.sort_by do |obj|
      atts.collect do |att|
        if down.include?(att)
          obj.send(att).rev
        else
          obj.send(att)
        end
      end
    end
  end

end
