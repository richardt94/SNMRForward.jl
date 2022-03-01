## set up McMC
using Distributed
nsamples, nchains, nchainsatone = 100001, 4, 1
Tmax = 2.50
addprocs(nchains)
##
@info "workers are $(workers())"
@everywhere using Distributed
@everywhere using transD_GP
@everywhere include("../SMRPI.jl")
## run McMC
@time transD_GP.main(opt, sounding, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
## close the worker pool
rmprocs(workers())
