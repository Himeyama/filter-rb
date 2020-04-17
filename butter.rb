#!/usr/bin/env ruby
require "narray"

class Filter
    def iirfilter(n, wn, rp=nil, rs=nil, btype="band", analog=false, ftype="butter", output="ba", fs=nil)
        btype, ftype, output = [btype, ftype, output].map(&:downcase)
        # p ftype
        wn = NArray[wn]
        # wn = asarray(wn)

        if fs
            raise "fs cannot be specified for an analog filter" if analog
            wn = 2 * wn / fs
        end

        # if fs is not None:
        #     if analog:
        #         raise ValueError("fs cannot be specified for an analog filter")
        #     Wn = 2*Wn/fs

        btype == "low" ? (btype = "lowpass") : raise("未対応")

        # begin
        #     btype = @@band_dict[btype]
        # rescue
        #     raise "'%s' is an invalid bandtype for filter." % btype
        # end

        ftype == "butter" ? (typefunc = "buttap") : raise("未対応")

        # try:
        #     btype = band_dict[btype]
        # except KeyError:
        #     raise ValueError("'%s' is an invalid bandtype for filter." % btype)

        # begin
        #     typefunc = @@filter_dict[ftype][0]
        # rescue
        #     raise "'%s' is not a valid basic IIR filter." % ftype
        # end

        # try:
        #     typefunc = filter_dict[ftype][0]
        # except KeyError:
        #     raise ValueError("'%s' is not a valid basic IIR filter." % ftype)

        # p ["ba", "zpk", "sos"].include?(output)
        raise "'%s' is not a valid output form." % output unless ["ba", "zpk", "sos"].include?(output)

        # if output not in ['ba', 'zpk', 'sos']:
        #     raise ValueError("'%s' is not a valid output form." % output)

        raise "passband ripple (rp) must be positive" if rp && rp < 0

        # if rp is not None and rp < 0:
        #     raise ValueError("passband ripple (rp) must be positive")

        raise "stopband attenuation (rs) must be positive" if rs && rs < 0

        # if rs is not None and rs < 0:
        #     raise ValueError("stopband attenuation (rs) must be positive")

        # Get analog lowpass prototype

        
        if typefunc == "buttap"
            z, p, k = buttap(n)
        else
            raise("未対応")
        end



        # if typefunc == buttap:
        #     z, p, k = typefunc(N)
        # elif typefunc == besselap:
        #     z, p, k = typefunc(N, norm=bessel_norms[ftype])
        # elif typefunc == cheb1ap:
        #     if rp is None:
        #         raise ValueError("passband ripple (rp) must be provided to "
        #                         "design a Chebyshev I filter.")
        #     z, p, k = typefunc(N, rp)
        # elif typefunc == cheb2ap:
        #     if rs is None:
        #         raise ValueError("stopband attenuation (rs) must be provided to "
        #                         "design an Chebyshev II filter.")
        #     z, p, k = typefunc(N, rs)
        # elif typefunc == ellipap:
        #     if rs is None or rp is None:
        #         raise ValueError("Both rp and rs must be provided to design an "
        #                         "elliptic filter.")
        #     z, p, k = typefunc(N, rp, rs)
        # else:
        #     raise NotImplementedError("'%s' not implemented in iirfilter." % ftype)

        # Pre-warp frequencies for digital filter design

        unless analog
            if numpy.any(Wn <= 0) or numpy.any(Wn >= 1)
                if fs
                    raise
                end
                raise
            end
            fs = 2.0
            warped = 2 * fs * tan(pi * Wn / fs)
        else
            warped = Wn
        end

        # if not analog:
        #     if numpy.any(Wn <= 0) or numpy.any(Wn >= 1):
        #         if fs is not None:
        #             raise ValueError("Digital filter critical frequencies "
        #                             "must be 0 < Wn < fs/2 (fs={} -> fs/2={})".format(fs, fs/2))
        #         raise ValueError("Digital filter critical frequencies "
        #                         "must be 0 < Wn < 1")
        #     fs = 2.0
        #     warped = 2 * fs * tan(pi * Wn / fs)
        # else:
        #     warped = Wn


        if ["lowpass", "highpass"].include?(btype)
            raise if wn.size != 1
            if btype == "lowpass"
                z, p, k = lp2lp_zpk(z, p, k, wo=warped)
            elsif btype == "highpass"
                z, p, k = lp2hp_zpk(z, p, k, wo=warped)
            end
        else
            raise "未対応"
        end

        # # transform to lowpass, bandpass, highpass, or bandstop
        # if btype in ('lowpass', 'highpass'):
        #     if numpy.size(Wn) != 1:
        #         raise ValueError('Must specify a single critical frequency Wn for lowpass or highpass filter')

        #     if btype == 'lowpass':
        #         z, p, k = lp2lp_zpk(z, p, k, wo=warped)
        #     elif btype == 'highpass':
        #         z, p, k = lp2hp_zpk(z, p, k, wo=warped)
        # elif btype in ('bandpass', 'bandstop'):
        #     try:
        #         bw = warped[1] - warped[0]
        #         wo = sqrt(warped[0] * warped[1])
        #     except IndexError:
        #         raise ValueError('Wn must specify start and stop frequencies for bandpass or bandstop filter')

        #     if btype == 'bandpass':
        #         z, p, k = lp2bp_zpk(z, p, k, wo=wo, bw=bw)
        #     elif btype == 'bandstop':
        #         z, p, k = lp2bs_zpk(z, p, k, wo=wo, bw=bw)
        # else:
        #     raise NotImplementedError("'%s' not implemented in iirfilter." % btype)

        z, p, k = bilinear_zpk(z, p, k, fs) unless analog

        # Find discrete equivalent if necessary
        # if not analog:
        #     z, p, k = bilinear_zpk(z, p, k, fs=fs)

        # Transform to proper out type (pole-zero, state-space, numer-denom)
        # if output == 'zpk':
        #     return z, p, k
        # elif output == 'ba':
        #     return zpk2tf(z, p, k)
        # elif output == 'sos':
        #     return zpk2sos(z, p, k)
        raise "未対応" unless output == "ba"
        [z, p, k]
    end

    def butter(n, wn, btype="low", analog=false, output="ba", fs=nil)
        iirfilter(n, wn, nil, nil, btype, analog, "butter", output, fs)
    end
end