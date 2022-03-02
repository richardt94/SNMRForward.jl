using transD_GP, PyPlot, MAT, DelimitedFiles
cd(@__DIR__)
include("../SMRPI.jl")
include("../ProcessingTools.jl")
using .SMRPI, .ProcessingTools
## Load the processed data from a GMR FID sounding
example_fid = matread("../../../example_data/FID_40ms.mat")
t = example_fid["time_fid"][:]
fid_qt = example_fid["coil_1_fid"]
##
V0, ϕ = get_sounding_curve(t, fid_qt)
## params for the sounding - some of these are stored in the MATLAB file
# but others (e.g. field inclination) are from site info for the survey
# or separate ASCII files
q = example_fid["pulse_moment"]
freq = example_fid["detect_frequency"]
L = 50
inclination = 43.9 * π/180 #degrees to radians
θ = 0 #loop oriented mag. north
Be = 2π*freq/γh
resist_data = readdlm("../../../example_data/1pm_res_profile.txt")
c = 1 ./ resist_data[2:end,1]
t = Vector{Float64}(resist_data[2:end-1,2])
σt = ConductivityModel(c,t)
## define a depth grid for the modelling and inversion
zstart = 1.
extendfrac, dz = 1.028, 1
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=60, showplot=true)
gcf()
##
F = MRSForward_square(L, zboundaries, q[:], inclination, 0, Be, σt)
sounding = newSMRSounding(V0[:], ϕ[:], F)
##