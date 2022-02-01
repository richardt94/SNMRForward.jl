module MyFFTUtils
using FFTW

function getf(fs, m)
    # NOT fftshifted, human readable
    if iseven(m)
        (-m/2:m/2-1)*fs/m
    else
        (-(m-1)/2:(m-1)/2)*fs/m
    end
end  

function BWfilt(f, f₀; n=2)
    # butterworth low pass (zero phase)
    L = 1 ./sqrt.(1 .+ (f/f₀).^(2n))
end 

function BWbandpassfilt(f, f1, f2; n=2)
    # two butterworth bandpass
    fmin, fmax = extrema([f1, f2])
    Linclude = BWfilt(f, fmax, n=n)
    Lexclude = BWfilt(f, fmin, n=n)
    Linclude - Lexclude
end

function expnotch(f, fc; n=2, λ=1)
    # exponential family notch
    L = exp.(-((f.-fc)/λ).^n)
end

function getquad(x)
    # apply the hilbert transform in frequency
    @assert all(isa.(x, Real))
    X = fft(x)
    m = length(X)
    f = ifftshift(getf(1, m))
    HX = -im*sign.(f).*X
    real(ifft(HX))
end

function unwrap(v, inplace=false)
    # currently assuming an array
    # could use a tolerance somewhere I guess
    unwrapped = inplace ? v : copy(v)
    for i in 2:length(v)
        while unwrapped[i] - unwrapped[i-1] >= pi
            unwrapped[i] -= 2pi
        end
        while unwrapped[i] - unwrapped[i-1] <= -pi
            unwrapped[i] += 2pi
        end
    end
    return unwrapped
end

export getf, BWFilt, BWbandpassfilt, expnotch, getquad, unwrap

end