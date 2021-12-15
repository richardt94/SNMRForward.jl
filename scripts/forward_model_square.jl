using SNMRForward, PyPlot
##
# square side length to match area of 50 m circular loop
L = sqrt(π*50^2/4)
log_zgrid = 0.5:0.05:2
zgrid = 10 .^ log_zgrid
## conductive half-space, at 2.5 kHz
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
pcolor(xgrid, zgrid, real.(slicex)', vmin = -0.01, vmax = 0.01)
gca().invert_yaxis()
sca(ax[2])
pcolor(xgrid, zgrid, real.(slicey)', vmin = -0.01, vmax = 0.01)
gca().invert_yaxis()
sca(ax[3])
pcolor(xgrid, zgrid, real.(slicez)', vmin = -0.01, vmax = 0.01)
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
pcolor(xgrid, zgrid, real.(slice_co)')
gca().invert_yaxis()
sca(ax[2])
pcolor(xgrid, zgrid, real.(slice_counter)')
gca().invert_yaxis()
display(gcf())
close(gcf())
##

q = 10
kz = SNMRForward.kernel_1d(q, ϕ, θ, ωl, Hx, Hy, Hz, xgrid)
##
figure()
plot(real.(kz),zgrid)
gca().invert_yaxis()
gca().set_yscale("log")
display(gcf())
##
μ = SNMRForward.μ0
θ = 0
Hparams = SNMRForward.co_counter_field.(Hx, Hy, Hz, ϕ, θ)

Hco = first.(Hparams)
ζ = last.(Hparams)
Hctr = reshape([a[2] for a in Hparams[:]], size(Hco)...)

q = 10
slice_co = Hco[midpoint, :, :]
αt = ((SNMRForward.γh * q * μ * slice_co) .% (2*π))/π * 180

figure()
cs = contourf(xgrid, zgrid, αt', levels=[0,45,90,135,210,225,270,315,360])
gca().invert_yaxis()
xlabel("distance from loop centre (m)")
ylabel("depth (m)")
colorbar(cs, label = "tipping angle (degrees)")
display(gcf())
# savefig("tipping.png")

##
qgrid = [0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
c = SNMRForward.ConductivityModel(σ, d)
F = SNMRForward.MRSForward_square(L, zgrid, qgrid, ϕ, θ, Be, c)
##
w = zeros(length(zgrid))
w[(zgrid .>= 10) .& (zgrid .<= 20)] .= 1
data = SNMRForward.forward(F,w)

figure()
plot(qgrid,real.(data))
display(gcf())
##

w = zeros(length(zgrid))
w[(zgrid .>= 10) .& (zgrid .<= 20)] .= 1
response = SNMRForward.forward(F,w)

fig, ax = subplots(1,3, figsize=(15,5))
sca(ax[1])
gca().ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plot(qgrid, real.(response))
title("Saturated layer 10 - 20 m")
ylabel("Response voltage (V)")
xlabel("Pulse moment (A s)")
display(gcf())
## 30 - 45 m
w = zeros(length(zgrid))
w[(zgrid .>= 30) .& (zgrid .<= 45)] .= 1

response = SNMRForward.forward(F,w)

sca(ax[2])
gca().ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plot(qgrid, real.(response))
title("Saturated layer 30 - 45 m")
xlabel("Pulse moment (A s)")
display(gcf())
## 60 - 80 m
w = zeros(length(zgrid))
w[(zgrid .>= 60) .& (zgrid .<= 80)] .= 1

response = SNMRForward.forward(F,w)

sca(ax[3])
gca().ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plot(qgrid, real.(response))
title("Saturated layer 60 - 80 m")
xlabel("Pulse moment (A s)")
display(gcf())


savefig("Forwards_square.png")
##
