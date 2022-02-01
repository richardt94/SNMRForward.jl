ts_txt = readlines("scripts/time_amplitude.txt")
##
t = parse.(Float64, first.(split.(ts_txt)))
ft = parse.(Float64, last.(split.(ts_txt)))
##
using PyPlot
##
figure()
plot(t,ft)
gcf()
##
using FFTW
##
ftspec = abs.(rfft(ft)).^2 / length(t)^2
freq = rfftfreq(length(t)) / (t[2]-t[1])
figure()
plot(freq, ftspec)
gca().set_yscale("log")
xlim([2090,2170])
# ylim([0,0.3])
xlabel("frequency")
ylabel("power")
gcf()
##
dominant_f = freq[ftspec .> 1e-4]
##
dominant_f[end-5:end]
##
using MAT
##
FID_data = matread("example_data/FID_40ms.mat")
##
for k in keys(FID_data)
    println(k)
end
##
any(FID_data["coil_1_fid"] .!= 0)
##
times = FID_data["time_fid"]
fidtimeseries = FID_data["coil_1_fid"]
##
times
##
ts_1 = fidtimeseries[:,1]
##
figure()
plot(times[:], ts_1)
xlim([0.015,0.1])
gcf()
##
close("all")
##
# FFT of the actual FID signal
dt = times[2] - times[1]
ndata = length(times)
spec = abs.(rfft(ts_1)).^2 / ndata^2
freq = rfftfreq(ndata) / dt
figure()
plot(freq, spec)
gca().set_yscale("log")
xlabel("frequency (Hz)")
ylabel("power (V^2)")
gcf()
##
# synthesise a signal
T = 0.15
ω = 2040
synth_timeseries = exp.(-times/T) .* sin.(2*pi*ω * times)
##
figure()
plot(times[:],synth_timeseries[:])
xlim([0.015,0.1])
gcf()
##

for qi = 1:size(fidtimeseries,2)
    if qi % 10 == 1
        fig, ax = subplots(10,2,figsize=(10,20), constrained_layout=true)
    end
    sca(ax[(qi-1)%10+1,1])
    tsi = fidtimeseries[:,qi]
    plot(times[:], tsi)
    spec = abs.(rfft(tsi)).^2 / ndata^2
    sca(ax[(qi-1)%10+1,2])
    plot(freq, spec)
    gca().set_yscale("log")
    xlabel("frequency (Hz)")
    ylabel("power (V^2)")
    if qi % 10 == 0
        savefig("FID_pulse_moments_$(qi ÷ 10).png")
    end
end
##
gcf()
##
close("all")