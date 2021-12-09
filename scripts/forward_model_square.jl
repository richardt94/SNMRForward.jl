using SNMRForward, PyPlot
##
# square side length to match area of 50 m circular loop
L = sqrt(π*50^2/4)
zgrid = 0:1:100.0
## conductive half-space, at 2 kHz
ωl = 2*π*2.5e3 #Hz, typical for Earth's field strength
Be = ωl/SNMRForward.γh
# Be = 0.000048
# ωl = SNMRForward.γh * 0.000048
d = Vector{Float64}()#[20.0,30.0]
σ = [0.001]#,0.1, 0.02]
##
# FFTs define position and momentum grids for us
# with a square loop (due to the properties of the DFT)
nxpoints = 128
Hx, Hy, Hz, xgrid, kgrid = SNMRForward.magfields_square(L, ωl, σ, d, 4*L, zgrid, nxpoints=nxpoints);
##
midpoint = nxpoints÷2 + 1
slicex = Hx[:,midpoint,:]
slicey = Hy[midpoint,:,:]
slicez = Hz[midpoint,:,:]
fig, ax = subplots(2,2, figsize=(10,10))
sca(ax[1])
pcolor(xgrid, zgrid, real.(slicex)', vmin=-3e-4, vmax=3e-4)
gca().invert_yaxis()
sca(ax[2])
pcolor(xgrid, zgrid, real.(slicey)', vmin=-3e-4, vmax=3e-4)
gca().invert_yaxis()
sca(ax[3])
pcolor(xgrid, zgrid, real.(slicez)', vmin=-3e-4, vmax=3e-4)
gca().invert_yaxis()

display(gcf())
close(gcf())
##
#compute co-rotating fields given Hx, Hy, Hz
θ = 0
ϕ = 13*π/36
Hparams = SNMRForward.co_counter_field.(Hx, Hy, Hz, ϕ, θ)

Hco = first.(Hparams)
ζ = last.(Hparams)

H_counter = reshape([a[2] for a in Hparams[:]], size(Hco)...)

##
slice_co = Hco[midpoint,:,:]
slice_counter = H_counter[midpoint,:,:]
fig, ax = subplots(1,2, figsize=(10,5))
sca(ax[1])
pcolor(xgrid, zgrid, real.(slice_co)', vmin=-3e-4, vmax=3e-4)
gca().invert_yaxis()
sca(ax[2])
pcolor(xgrid, zgrid, real.(slice_counter)', vmin=-3e-4, vmax=3e-4)
gca().invert_yaxis()
display(gcf())
close(gcf())
##
