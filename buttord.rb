#!/usr/bin/env ruby
require "narray"

class Filter
    def buttord(wp, ws, gpass, gstop, analog=false, fs=nil)
        wp = NArray[wp]
        ws = NArray[ws]
        raise "fs cannot be specified for an analog filter" if fs && analog
        wp = 2 * wp / fs if fs
        ws = 2 * ws / fs if fs
        passb = NMath::tan(Math::PI * wp / 2.0) unless analog
        stopb = NMath::tan(Math::PI * ws / 2.0) unless analog
        passb = wp * 1.0 if analog
        stopb = ws * 1.0 if analog
        nat = stopb / passb
        nat = nat.abs.min
        gSTOP = 10 ** (0.1 * gstop.abs)
        gPASS = 10 ** (0.1 * gpass.abs)
        ord = (NMath::log10((gSTOP - 1.0) / (gPASS - 1.0)) / (2 * NMath::log10(nat))).ceil.to_i
        w0 = (gPASS - 1.0) ** (-1.0 / (2.0 * ord))
        wN = w0 * passb
        wn = analog ? wN : (2.0 / Math::PI) * NMath::atan(wN)
        wn = wn[0] if wn.size == 1
        wn = wn * fs / 2 if fs
        [ord, wn]
    end
end