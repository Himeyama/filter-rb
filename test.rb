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

class Filter
    def buttord(wp, ws, gpass, gstop, analog=false, fs=nil)
        wp = NArray[wp]
        ws = NArray[ws]
        if fs
            raise "fs cannot be specified for an analog filter" if analog
            wp = 2 * wp / fs
            ws = 2 * ws / fs
        end
        filter_type = 2 * (wp.size - 1) + 1
        filter_type += 1 if wp[0] >= ws[0]
        unless analog
            passb = NMath::tan(Math::PI * wp / 2.0)
            stopb = NMath::tan(Math::PI * ws / 2.0)
        else
            passb = wp * 1.0
            stopb = ws * 1.0
        end
        if filter_type == 1
            nat = stopb / passb
        elsif filter_type == 2
            nat = passb / stopb
        elsif filter_type == 3
            wp0 = optimize.fminbound(
                band_stop_obj, 
                passb[0], 
                stopb[0] - 1e-12,
                args=[0, passb, stopb, gpass, gstop, 'butter'], 
                disp=0
            )
            passb[0] = wp0
            wp1 = optimize.fminbound(
                band_stop_obj, 
                stopb[1] + 1e-12, 
                passb[1],
                args=[1, passb, stopb, gpass, gstop, 'butter'],
                disp=0
            )
            passb[1] = wp1
            nat = ((stopb * (passb[0] - passb[1])) / (stopb ** 2 - passb[0] * passb[1]))
        elsif filter_type == 4
            nat = ((stopb ** 2 - passb[0] * passb[1]) / (stopb * (passb[0] - passb[1])))
        end
        nat = nat.abs.min
        gSTOP = 10 ** (0.1 * gstop.abs)
        gPASS = 10 ** (0.1 * gpass.abs)
        ord = (NMath::log10((gSTOP - 1.0) / (gPASS - 1.0)) / (2 * NMath::log10(nat))).ceil.to_i
        w0 = (gPASS - 1.0) ** (-1.0 / (2.0 * ord))
        if filter_type == 1
            wN = w0 * passb
        elsif filter_type == 2
            wN = passb / w0
        elsif filter_type == 3
            wN = numpy.zeros(2, float)
            discr = sqrt((passb[1] - passb[0]) ** 2 +
                        4 * w0 ** 2 * passb[0] * passb[1])
            wN[0] = ((passb[1] - passb[0]) + discr) / (2 * w0)
            wN[1] = ((passb[1] - passb[0]) - discr) / (2 * w0)
            wN = numpy.sort(abs(wN))
        elsif filter_type == 4
            w0 = numpy.array([-w0, w0], float)
            wN = (-w0 * (passb[1] - passb[0]) / 2.0 +
                sqrt(w0 ** 2 / 4.0 * (passb[1] - passb[0]) ** 2 +
                    passb[0] * passb[1]))
            wN = numpy.sort(abs(wN))
        else
            raise "Bad type: %s" % filter_type
        end
        wn = analog ? wN : (2.0 / Math::PI) * NMath::atan(wN)
        wn = wn[0] if wn.size == 1
        wn = wn * fs / 2 if fs
        [ord, wn]
    end
end