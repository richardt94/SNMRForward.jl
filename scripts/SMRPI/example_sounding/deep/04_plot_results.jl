using Statistics
##

close("all")
transD_GP.getchi2forall(opt)
ax = gcf().axes;
# χ² = length(sounding.V0)
# ax[2].plot(xlim(), [χ², χ²], "--", color="gray")
# ax[2].set_ylim(χ²-10, χ²+10)
gcf()
##
istothepow = false
@assert !(linearsat & istothepow)
opt.xall[:] .= zboundaries
transD_GP.plot_posterior(sounding, opt, burninfrac=0.5, figsize=(10,6), qp1=0.2, qp2=0.8, nbins=50, istothepow=istothepow, cmappdf="winter",
    vmaxpc=1.0, pdfnormalize=false, plotmean=false, lwidth=1)
ax = gcf().axes
linearsat ? ax[1].set_xlabel("fractional water content") : ax[1].set_xlabel("log\$_{10}\$ water content") 
linearsat || ax[1].plot(istothepow ? wc_gmr : log10.(wc_gmr), z_gmr, "w-")
linearsat && ax[1].plot(wc_gmr, z_gmr, "w-")
gcf()
##
# nuisance histograms
transD_GP.plot_posterior(sounding, optn, burninfrac=0.5, nbins=50)
gcf()

##
F = transD_GP.assembleTat1(opt, :U, temperaturenum=1)
est_σ2 = exp.(1/length(sounding.V0) * F)/(2*length(sounding.V0))
est_σ = sqrt.(est_σ2)

mean(est_σ)