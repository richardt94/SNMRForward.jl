include("MyFFTUtils.jl")
using .MyFFTUtils, PyPlot, FFTW, Random, Statistics
##
Random.seed!(10)
fs = 10_000 # sampling frequency
t = 0:1/fs:0.5 # time vector
f = 2_000 # Larmor frequency
E0, df, ϕ = 10., 2., deg2rad(160)
T2 = 0.5 # Decay time
x = E0*cos.(2*pi*(f+df)*t .+ ϕ).*exp.(-t/T2) # Synthetic NMR signal
freqs = getf(fs, length(t)) # from MyFFTUtils get the frequency vector for FFT
X = fft(x)
## create a noisy synthetic and analytic signal with a dead time
tdead = 0.005 
noisefrac = 0.3 # fraction of E0
idx = findfirst(t.>=tdead)
xnoisy = x[idx+1:end] + noisefrac*E0*randn(length(x)-idx)
ynoisy = getquad(xnoisy) # Hilbert transform in MyFFTUtils
znoisy = xnoisy + im*ynoisy # Form analytic signal
tnoisy = t[idx+1:end]
XN = fft(xnoisy)
freqnoisy = getf(fs, length(t)-idx)
## plot noisy and true values
fig = figure()
s1 = subplot(211)
s1.plot(t, x)
s1.plot(t[idx+1:end], xnoisy, alpha=0.6)
xlabel("Time s")
s2 = subplot(212)
s2.semilogy(freqs, abs.(fftshift(X)))
s2.semilogy(freqnoisy, abs.(fftshift(XN)), alpha=0.2, "-o")
xlabel("Frequency Hz")
fig.tight_layout()
gcf()
## now get initial estimates: initial value and decay T2
Gₐ = [-tnoisy ones(size(tnoisy))]
dₐ = log.(abs.(znoisy)) 
mₐ = (Gₐ'Gₐ)\(Gₐ'dₐ)
@show T2est = 1/mₐ[1]
@show E0est = exp(mₐ[2])
## next df and ϕ
dᵩ = unwrap(angle.(znoisy)) - 2pi*f*tnoisy
Gᵩ = [2*pi*tnoisy ones(length(tnoisy))]
mᵩ = (Gᵩ'Gᵩ)\(Gᵩ'dᵩ)
@show dfest = mᵩ[1]
@show ϕest = mod(rad2deg(mᵩ[2]), 360)
## plot estimates on new figure
xest = E0est*cos.(2*pi*(f+dfest)*t .+ deg2rad(ϕest)).*exp.(-t/T2est)
Xest = fft(xest)
fig = figure()
s1 = subplot(211)
s1.plot(t[idx+1:end], xnoisy, color="orange")
s1.plot(t, x, label="original")
s1.plot(t, xest, "-r", label="estimated")
legend()
xlabel("Time s")
s2 = subplot(212)
s2.semilogy(freqnoisy, abs.(fftshift(XN)), "-o", color="orange")
s2.semilogy(freqs, abs.(fftshift(X)))
s2.semilogy(freqs, abs.(fftshift(Xest)), "-r")
xlabel("Frequency Hz")
fig.tight_layout()
gcf()
## print out results
@info "True values"
@info "df: $df Hz ϕ: $(rad2deg(ϕ)) deg E0: $E0 T2: $T2 sec"
@info "Estimated"
@info "df: $dfest Hz ϕ: $ϕest deg E0: $E0est T2: $T2est sec"