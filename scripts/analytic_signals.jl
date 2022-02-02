using PyPlot, FFTW, Random

##
c = 10 * rand()
ϕ = 2*pi*rand() - pi
ω = 1
##
t = 16*pi*(0:0.001:1)
sig_real = c*cos.(ω*t .+ ϕ)
##
figure()
plot(t,sig_real)
gcf()
##
sigf = fft(sig_real)
freq = fftfreq(length(t))
anaf = ((freq .>= 0) .+ (freq .> 0)) .* sigf
anat = ifft(anaf)
##
figure()
plot(t,sig_real,t,real.(anat))
gcf()
##
figure()
plot(t,angle.(anat), t, ((ω*t .+ ϕ) .+ π) .% (2 * pi) .- pi)
gcf()
##
using MAT
##
example_fid = matread("example_data/FID_40ms.mat")
t = example_fid["time_fid"][:]
fid_qt = example_fid["coil_1_fid"]
##
figure()
pcolor(fid_qt)
gcf()
##
fid_ts = fid_qt[:,1]
figure()
plot(fid_ts)
gcf()
##
function analytic_signal(signal)
    sigf = fft(signal)
    freq = fftfreq(length(signal))
    anaf = ((freq .>= 0) .+ (freq .> 0)) .* sigf
    ifft(anaf)
end
##
figure()
plot(t, abs.(analytic_signal(fid_ts)))
gcf()
##
log_amp = log.(abs.(analytic_signal(fid_ts)))
##
figure()
plot(t, log_amp)
gcf()
##
using LinearAlgebra
##
G = hcat(-t, ones(length(t)))
invT2, logE0 = pinv(G)*log_amp
##
T2 = 1/invT2
E0 = exp(logE0)
##
model = logE0 .- invT2*t
##
figure()
plot(t,model,t,log_amp)
gcf()
##
model_linear = E0 * exp.(-t/T2)
figure()
plot(t,model_linear,t,abs.(analytic_signal(fid_ts)))
gcf()
## first-order error propagation
using Statistics
##
amp = abs.(analytic_signal(fid_ts))
vard = var(amp .- model_linear)
σd = sqrt(vard)
##
σlog = σd ./ amp
Cd = diagm(σlog.^2)
Cm = pinv(G)*Cd*transpose(pinv(G))
σm = sqrt.(diag(Cm))
σT2 = T2^2 * σm[1]
σE0 = E0*σm[2]
@info "T2 = $T2 ± $σT2"
@info "E0 = $E0 ± $σE0"
##
function est_E0_T2(t, timeseries)
    amp = abs.(analytic_signal(timeseries))
    log_amp = log.(amp)

    G = hcat(-t, ones(length(t)))
    regr = pinv(G)
    est = regr*log_amp
    T2 = 1/est[1]
    E0 = exp(est[2])
    modelled_amp = E0 * exp.(-t/T2)
    vard = var(amp .- modelled_amp)
    σlog = sqrt(vard) ./ amp
    Cd = diagm(σlog.^2)
    Cm = regr * Cd * regr'
    σm = sqrt.(diag(Cm))
    σT2 = T2^2 * σm[1]
    σE0 = E0*σm[2]

    (E0, T2, σE0, σT2)
end
##
est_E0_T2(t, fid_ts)
##
estimates = mapslices(x -> est_E0_T2(t,x), fid_qt; dims=1)   
##
sounding = first.(estimates)
sounding_err = [e[3] for e in estimates]
##
q = example_fid["pulse_moment"]
##
figure()
xlabel("pulse moment (A s)")
ylabel("E0 (V)")
errorbar(q[:],sounding[:], yerr=2*sounding_err, fmt=".")
gcf()
##
figure()
xlabel("pulse moment (A s)")
ylabel("E0 (V)")
scatter(q[:], 2*sounding_err)
gcf()
##
T2 = [e[2] for e in estimates]
ω = 2*π*2040
ϕ = rand(1,20) * 2 * π
##
noise_mag = 100e-9
synth_data = sounding .* cos.(ω*t .+ ϕ) .* exp.(-t./T2)
noisy_synth = noise_mag .* randn(size(synth_data)) .+ synth_data
##
estimates = mapslices(x -> est_E0_T2(t,x), noisy_synth; dims=1)   
##
sounding2 = first.(estimates)
sounding2_err = [e[3] for e in estimates]
figure()
xlabel("pulse moment (A s)")
ylabel("E0 (V)")
errorbar(q[:],sounding2[:], yerr=2*sounding2_err, fmt=".")
gcf()
## estimate frequency and phase
function unwrap!(phase)
    for i in 2:length(phase)
        while phase[i] - phase[i-1] >= pi
            phase[i] -= 2pi
        end
        while phase[i] - phase[i-1] <= -pi
            phase[i] += 2pi
        end
    end
end
##
arg = angle.(analytic_signal(fid_ts))
unwrap!(arg)
##
figure()
plot(t,arg)
gcf()
##
function est_ω_ϕ(t, timeseries)
    arg = angle.(analytic_signal(timeseries))
    unwrap!(arg)
    G = hcat(t, ones(length(t)))
    regr = pinv(G)
    ω, ϕ = regr*arg
    vard = var(ω*t .+ ϕ .- arg)
    Cm = vard * regr * regr'
    σω, σϕ = sqrt.(diag(Cm))

    (ω, ϕ, σω, σϕ)
end

##
E0, T2, _, _ = est_E0_T2(t, fid_ts)
ω, ϕ, _, _ = est_ω_ϕ(t, fid_ts)
##
modelled_signal = E0 * cos.(ω*t .+ ϕ) .* exp.(-t/T2)
##
figure()
plot(t, fid_ts, t, modelled_signal)
gcf()
##
function full_estimates(t, ts)
    (est_E0_T2(t,ts)..., est_ω_ϕ(t,ts)...)
end
##
estimates = mapslices(x -> full_estimates(t,x), fid_qt; dims=1)
##
E0 = [e[1] for e in estimates]
ω = [e[5] for e in estimates]
ϕ = [e[6] for e in estimates]
##
figure()
plot(q[:], ϕ[:] .% (2π))
gcf()
##
V0 = E0 .* exp.(ϕ*im)
figure()
plot(q[:], abs.(V0[:]))
gcf()
##
using SNMRForward
##
# read GMR inversion result
using DelimitedFiles
##
GMR_res = readdlm("example_data/conductive_earth_inversion_FID_40ms/conductive_earth_inversion_FID_40ms_1d_inversion.txt")
##
figure()
pcolor(GMR_res, vmin=0, vmax=1)
gcf()
##
size(GMR_res)
##
# wc_gmr = sum(GMR_res[1:end-1,12:end], dims=2)
wc_gmr = GMR_res[1:end-1,4]
z_gmr = (GMR_res[1:end-1, 1] .+ GMR_res[1:end-1, 1])/2
q_gmr = example_fid["pulse_moment"][:]
ω_gmr = 2π*sum(GMR_res[1:end-1, 6]) / (size(GMR_res, 1) - 1)
## just go with a resistive earth for now
σ = [0.004]
d = Vector{Float64}()
σd = SNMRForward.ConductivityModel(σ, d)
## Guess a 100 m square loop
ϕ = -43.9*π/180
Be = ω_gmr / SNMRForward.γh
F = SNMRForward.MRSForward_square(50, z_gmr, q_gmr, ϕ, 0, Be, σd)

##
Vq = SNMRForward.forward(F, wc_gmr[:])
figure()
plot(q_gmr, real.(Vq), q_gmr, abs.(V0[:]))
legend(["forward of inverted WC", "V0 fitted data"])
gcf()

## load the site resistivity profile
resist_data = readdlm("example_data/1pm_res_profile.txt")
c = 1 ./ resist_data[2:end,1]
d = Vector{Float64}(resist_data[2:end-1,2])
##
σd = SNMRForward.ConductivityModel(c,d)

F = SNMRForward.MRSForward_square(50, z_gmr, q_gmr, ϕ, 0, Be, σd)
##
Vq = SNMRForward.forward(F, wc_gmr[:])
figure()
plot(q_gmr, abs.(Vq), q_gmr, abs.(V0[:]))
legend(["forward of inverted WC", "V0 fitted data"])
gcf()
##