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
        # _validate_gpass_gstop(gpass, gstop)

        # wp = atleast_1d(wp)
        # ws = atleast_1d(ws)

        wp = NArray[wp]
        ws = NArray[ws]

        if fs
            if analog
                raise "fs cannot be specified for an analog filter"
            end
            wp = 2 * wp / fs
            ws = 2 * ws / fs
        end

        # if fs is not None:
        #     if analog:
        #         raise ValueError("fs cannot be specified for an analog filter")
        #     wp = 2*wp/fs
        #     ws = 2*ws/fs

        filter_type = 2 * (wp.size - 1) + 1

        # filter_type = 2 * (len(wp) - 1)
        # filter_type += 1

        if wp[0] >= ws[0]
            filter_type += 1
        end

        # if wp[0] >= ws[0]:
        #     filter_type += 1

        # Pre-warp frequencies for digital filter design

        unless analog
            passb = NMath::tan(Math::PI * wp / 2.0)
            stopb = NMath::tan(Math::PI * ws / 2.0)
        else
            passb = wp * 1.0
            stopb = ws * 1.0
        end

        # if not analog:
        #     passb = tan(pi * wp / 2.0)
        #     stopb = tan(pi * ws / 2.0)
        # else:
        #     passb = wp * 1.0
        #     stopb = ws * 1.0

        # nat = nil
        # p filter_type, passb, stopb
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

        # if filter_type == 1:            # low
        #     nat = stopb / passb
        # elif filter_type == 2:          # high
        #     nat = passb / stopb
        # elif filter_type == 3:          # stop
        #     wp0 = optimize.fminbound(band_stop_obj, passb[0], stopb[0] - 1e-12,
        #                             args=(0, passb, stopb, gpass, gstop,
        #                                 'butter'),
        #                             disp=0)
        #     passb[0] = wp0
        #     wp1 = optimize.fminbound(band_stop_obj, stopb[1] + 1e-12, passb[1],
        #                             args=(1, passb, stopb, gpass, gstop,
        #                                 'butter'),
        #                             disp=0)
        #     passb[1] = wp1
        #     nat = ((stopb * (passb[0] - passb[1])) / (stopb ** 2 - passb[0] * passb[1]))
        # elif filter_type == 4:          # pass
        #     nat = ((stopb ** 2 - passb[0] * passb[1]) / (stopb * (passb[0] - passb[1])))

        nat = nat.abs.min

        # nat = min(abs(nat))

        gSTOP = 10 ** (0.1 * gstop.abs)
        gPASS = 10 ** (0.1 * gpass.abs)
        ord = (NMath::log10((gSTOP - 1.0) / (gPASS - 1.0)) / (2 * NMath::log10(nat))).ceil.to_i

        # GSTOP = 10 ** (0.1 * abs(gstop))
        # GPASS = 10 ** (0.1 * abs(gpass))
        # ord = int(ceil(log10((GSTOP - 1.0) / (GPASS - 1.0)) / (2 * log10(nat))))

        # Find the Butterworth natural frequency WN (or the "3dB" frequency")
        # to give exactly gpass at passb.

        w0 = (gPASS - 1.0) ** (-1.0 / (2.0 * ord))

        # 略
        # try:
        #     W0 = (GPASS - 1.0) ** (-1.0 / (2.0 * ord))
        # except ZeroDivisionError:
        #     W0 = 1.0
        #     print("Warning, order is zero...check input parameters.")

        # now convert this frequency back from lowpass prototype
        # to the original analog filter


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



        # if filter_type == 1:  # low
        #     WN = W0 * passb
        # elif filter_type == 2:  # high
        #     WN = passb / W0
        # elif filter_type == 3:  # stop
        #     WN = numpy.zeros(2, float)
        #     discr = sqrt((passb[1] - passb[0]) ** 2 +
        #                 4 * W0 ** 2 * passb[0] * passb[1])
        #     WN[0] = ((passb[1] - passb[0]) + discr) / (2 * W0)
        #     WN[1] = ((passb[1] - passb[0]) - discr) / (2 * W0)
        #     WN = numpy.sort(abs(WN))
        # elif filter_type == 4:  # pass
        #     W0 = numpy.array([-W0, W0], float)
        #     WN = (-W0 * (passb[1] - passb[0]) / 2.0 +
        #         sqrt(W0 ** 2 / 4.0 * (passb[1] - passb[0]) ** 2 +
        #             passb[0] * passb[1]))
        #     WN = numpy.sort(abs(WN))
        # else:
        #     raise ValueError("Bad type: %s" % filter_type)

        unless analog
            wn = (2.0 / Math::PI) * NMath::atan(wN)
        else
            wn = wN
        end

        # if not analog:
        #     wn = (2.0 / pi) * arctan(WN)
        # else:
        #     wn = WN

        if wn.size == 1
            wn = wn[0]
        end

        # if len(wn) == 1:
        #     wn = wn[0]

        if fs
            wn = wn * fs / 2
        end

        # if fs is not None:
        #     wn = wn*fs/2

        [ord, wn]
    end
end