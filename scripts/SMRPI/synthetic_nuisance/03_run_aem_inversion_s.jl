## set up McMC
using Distributed
nsamples, nchains, nchainsatone = 100001, 4, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere using Distributed
@everywhere using transD_GP
@everywhere include("../SMRPI.jl")
## run McMC
@time begin
    if sounding.amponly
        transD_GP.main(opt, sounding, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
    else    
        transD_GP.main(opt, optn, sounding, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
    end
end    
## close the worker pool
rmprocs(workers())
