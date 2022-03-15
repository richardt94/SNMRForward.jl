using Statistics
##
close("all")
transD_GP.getchi2forall(opt)
if !noise_mle
    ax = gcf().axes;
    χ² = amponly ? length(sounding.V0)/2 : length(sounding.V0)
    ax[2].plot(xlim(), [χ², χ²], "--", color="gray")
    ax[2].set_ylim(χ²-10, χ²+10)
end
##
istothepow = false
@assert !(linearsat & istothepow)
opt.xall[:] .= zboundaries
transD_GP.plot_posterior(sounding, opt, burninfrac=0.5, figsize=(10,6), qp1=0.05, qp2=0.95, nbins=50, istothepow=istothepow, cmappdf="winter",
    vmaxpc=1.0, pdfnormalize=false, plotmean=false, lwidth=1)
ax = gcf().axes
linearsat ? ax[1].set_xlabel("fractional water content") : ax[1].set_xlabel("log\$_{10}\$ water content") 
linearsat || ax[1].step(istothepow ? w : log10.(w), zboundaries, "w-")
linearsat && ax[1].step(w, zboundaries, "w-")
## nuisance histograms
transD_GP.plot_posterior(sounding, optn, burninfrac=0.5, nbins=50, figsize=(5,4))
## swarm plots
if amponly
    SMRPI.plot_model_field(sounding, opt, decfactor=10, lcolor="k", modelalpha=0.08)
else    
    SMRPI.plot_model_field(sounding, opt, optn, decfactor=10, lcolor="k", modelalpha=0.08)
end    
##
if noise_mle
    ndata = amponly ? length(sounding.V0) : 2*length(sounding.V0)
    F = transD_GP.assembleTat1(opt, :U, temperaturenum=1)
    est_σ2 = exp.(2/ndata * F)/ndata
    est_σ = sqrt.(est_σ2)
    @info mean(est_σ)
end