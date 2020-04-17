#!/usr/bin/env ruby
require "narray"

module NMath
    def sinc(x)
        sinc = NMath::sin(Math::PI * x) / (Math::PI * x)
        x.where2[1].each do |i|
            sinc[i] = 1.0
        end
        sinc
    end

    module_function :sinc
end
