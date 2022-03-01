using Statistics
##

close("all")
transD_GP.getchi2forall(opt)
ax = gcf().axes;
χ² = length(sounding.V0)
ax[2].plot(xlim(), [χ², χ²], "--", color="gray")
ax[2].set_ylim(χ²-10, χ²+10)
gcf()
##
opt.xall[:] .= zboundaries
transD_GP.plot_posterior(sounding, opt, burninfrac=0.5, figsize=(10,6), qp1=0.05, qp2=0.95, nbins=50, vmaxpc=1.0)
ax = gcf().axes
ax[1].step(log10.(w[2:end]), zboundaries[2:end], color="w", linewidth=2)
ax[1].step(log10.(w[2:end]), zboundaries[2:end], color="g", linestyle="--")

ax[1].set_xlabel("log\$_{10}\$ water content")

gcf()
##
# nuisance histograms
transD_GP.plot_posterior(sounding, optn, burninfrac=0.5, nbins=50)
gcf()

##
F = transD_GP.assembleTat1(opt, :U, temperaturenum=1)
est_σ2 = exp.(2/38 * F)/38
est_σ = sqrt.(est_σ2)

mean(est_σ)