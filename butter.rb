#!/usr/bin/env ruby
require 'narray'

class Filter
    def _relative_degree(z, p)
        degree = p.size - z.size
        raise if degree < 0
        degree
    end

    def bilinear_zpk(z, p, k, fs)
        z = NArray[z]
        p = NArray[p]
        degree = _relative_degree(z, p)
        fs2 = 2.0 * fs
        z_z = (fs2 + z) / (fs2 - z)
        p_z = (fs2 + p) / (fs2 - p)
        z_z = NArray[*z_z.to_a.append(*(-(NArray.float(degree) + 1)).to_a)]
        fs2dz = (fs2 - z).prod == NArray[] ? 1.0 : (fs2 - z).prod
        k_z = k * (fs2dz / (fs2 - p).prod).real
        [z_z, p_z, k_z]
    end

    def lp2lp_zpk(z, p, k, wo=1.0)
        z = NArray[z]
        p = NArray[p]
        wo = wo.to_f
        degree = _relative_degree(z, p)
        z_lp = wo * z
        p_lp = wo * p
        k_lp = k * wo ** degree
        [z_lp, p_lp, k_lp]
    end

    def buttap(n)
        raise "Filter order must be a nonnegative integer" unless n.to_i.abs == n
        z = NArray[]
        m = NArray.sfloat(n).indgen!(1 - n, 2)
        p = - NMath.exp(NArray[*Array.new(m.size){|i| 1i}] * Math::PI * m / (2 * n))
        [z, p, 1]
    end

    def iirfilter(n, wn, rp=nil, rs=nil, btype="band", analog=false, ftype="butter", output="ba", fs=nil)
        btype, ftype, output = [btype, ftype, output].map(&:downcase)
        wn = NArray[wn]
        if fs
            raise "fs cannot be specified for an analog filter" if analog
            wn = 2 * wn / fs
        end
        btype == "low" ? (btype = "lowpass") : raise("未対応")
        ftype == "butter" ? (typefunc = "buttap") : raise("未対応")
        raise "'%s' is not a valid output form." % output unless ["ba", "zpk", "sos"].include?(output)
        raise "passband ripple (rp) must be positive" if rp && rp < 0
        raise "stopband attenuation (rs) must be positive" if rs && rs < 0       
        if typefunc == "buttap"
            z, p, k = buttap(n)
        else
            raise("未対応")
        end
        unless analog
            raise if wn.le(0).where.size != 0 || wn.ge(1).where.size != 0
            fs = 2.0
            warped = 2 * fs * NMath::tan(Math::PI * wn / fs)
        else
            warped = wn
        end
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
        z, p, k = bilinear_zpk(z, p, k, fs) unless analog
        raise "未対応" unless output == "ba"
        [z, p, k]
    end

    def butter(n, wn, btype="low", analog=false, output="ba", fs=nil)
        iirfilter(n, wn, nil, nil, btype, analog, "butter", output, fs)
    end
end