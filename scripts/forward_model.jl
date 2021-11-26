using Revise, SNMRForward, PyPlot

R = 50
rgrid = 0.1:0.3:4*R
zgrid = 0.1*R:1:100

## conductive half-space, at 2 kHz
ωl = 2*π*2.5e3 #Hz, typical for Earth's field strength
Be = ωl/SNMRForward.γh
# Be = 0.000048
# ωl = SNMRForward.γh * 0.000048
d = [10.0,20.0]
σ = [0.001,0.05, 0.02]

# #define k grid for j0 and j1 kernels

# #j0        
# kj0 = reduce(hcat, SNMRForward.Filter_base_801/r for r in rgrid)
# rj0 = rgrid' .* ones(length(SNMRForward.Filter_base_801), length(rgrid))
# #j1
# kj1 = SNMRForward.Filter_base_801 / R

# ## determine whether to use j0 or j1 for hz at each r
# use_j1 = [r < R for r in rgrid]
# j1_inds = cumsum(use_j1)

# ## compute layer responses to propagate fields
# Bj0, αj0 = SNMRForward.responses(σ, d, kj0[:], ωl)
# Bj1, αj1 = SNMRForward.responses(σ, d, kj1, ωl)

# ## propagation coefficients for phi and phi' (Hz and Hr)
# phicj0, phipj0 = SNMRForward.phi_coeffs(Bj0, αj0, d, zgrid)
# phicj1, phipj1 = SNMRForward.phi_coeffs(Bj1, αj1, d, zgrid)

# ## calculate z = 0 for j0 and j1 kernels

# phif_j0 = SNMRForward.phi_free.(kj0[:], 0, R)

# phi0_j0 = 2 * kj0[:]./(kj0[:] .+ Bj0) .* phif_j0

# ##
# rj1 = rgrid[use_j1]
# ##
# phif_j1 = reduce(hcat, SNMRForward.phi_free_j1k.(kj1, 0, R, r) for r in rj1)

# phi0_j1 = 2 * kj1 ./ (kj1 .+ Bj1) .* phif_j1

# ## propagate through halfspace
# phiz_j0 = phi0_j0 .* phicj0
# phipz_j0 = phi0_j0 .* phipj0

# ##
# phiz_j1 = cat([phi0_j1[:,ir] .* phicj1 for ir=1:size(phi0_j1,2)]...; dims=3)

# ##low-pass filters
# filterspacing = SNMRForward.Filter_base_801[2] / SNMRForward.Filter_base_801[1]
# kcut(r) = pi/((filterspacing-1)*r)
# lowpass_k(k,r) = 1/(1 + (k/kcut(r))^3)

# # Hz_kernel_j1 =  lowpass_k.(kj1,R) .* kj1 .^ 3 .* permutedims(phiz_j1, (1,3,2))
# Hz_kernel_j1 = kj1 .^ 3 .* permutedims(phiz_j1, (1,3,2))


# ## compute Hz and Hr field kernels for each k value
# Hz_kernel = reshape(kj0[:].^3 .* phiz_j0, size(kj0)..., length(zgrid))
# # note the Hr kernel is actually in j1 space
# # due to the hankel transform property used to derive it
# Hr_kernel = reshape(-kj0[:].^2 .* phipz_j0, size(kj0)..., length(zgrid))

# # recover z and r fields on z and r grid
# Hz = zeros(ComplexF64, length(rgrid), length(zgrid))
# Hr = zeros(ComplexF64, length(rgrid), length(zgrid))
# for ir in 1:length(rgrid), iz in 1:length(zgrid)
#     if use_j1[ir]
#         ir1 = j1_inds[ir]
#         Hz[ir,iz] = 1/R * SNMRForward.Filter_J1_801' * Hz_kernel_j1[:,ir1,iz]
#     else
#         Hz[ir, iz] = 1/rgrid[ir] * SNMRForward.Filter_J0_801' * Hz_kernel[:,ir,iz]
#     end
#     Hr[ir, iz] = 1/rgrid[ir] * SNMRForward.Filter_J1_801' * Hr_kernel[:,ir,iz]
# end

##
(Hz,Hr) = SNMRForward.magfields(R,ωl,σ,d,rgrid,zgrid)

## plot the real part of the fields
fig, ax = subplots(1,2)
sca(ax[1])
pcolor(rgrid, zgrid, real.(Hz)', vmin = -0.01, vmax=0.01)
gca().invert_yaxis()
sca(ax[2])
pcolor(rgrid, zgrid, real.(Hr)', vmin = -0.01, vmax=0.01)
gca().invert_yaxis()
display(gcf())
close(gcf())
## say the mag. field inclination is 60 degrees. calculate corotating part and phase lag
ϕ = 12*π/36
Hfield_params = SNMRForward.co_counter_field.(Hz, Hr, ϕ, -π/2)

Hco = first.(Hfield_params)
ζ = last.(Hfield_params)


H_counter = reshape([a[2] for a in Hfield_params[:]], size(Hco)...)

##
fig, ax = subplots(1,2,figsize=(10,8))
sca(ax[1])
ylabel("z (m)")
xlabel("r (m)")
title("co-rotating field")
contourf(rgrid, zgrid, log.(Hco)', [-12,-7,-6.5,-6.0,-5.5,-5.0,-4.0,-3.0,-1.0])
gca().invert_yaxis()
sca(ax[2])
title("counter-rotating field")
xlabel("r (m)")
cs = contourf(rgrid, zgrid, log.(H_counter)', [-12,-7,-6.5,-6.0,-5.5,-5.0,-4.0,-3.0,-1.0])
gca().invert_yaxis()

colorbar(cs,location="bottom", ax=ax, label = "log H")

display(gcf())
savefig("co_counter_compare.png")
close(gcf())

##

figure()
cs = contourf(rgrid, zgrid, ζ')
gca().invert_yaxis()
colorbar(cs)
display(gcf())

## plot up the co- and counter-rotating B fields
# H -> B in tesla or gauss
# for a 300 A current, normalised to match figure 1 of weichman
μ_G = 4*π * 1e-3
norm_factor = 300 * μ_G / 0.299895
fig, ax = subplots(1,2,figsize=(10,8))
sca(ax[1])
ylabel("z (m)")
xlabel("r (m)")
title("co-rotating field")
contourf(rgrid, zgrid, log.(norm_factor * Hco)', [-11,-6,-5.5,-5.0,-4.5,-4.0,-3.0,-2.0,-0.0], cmap="jet")
gca().invert_yaxis()
sca(ax[2])
title("counter-rotating field")
xlabel("r (m)")
cs = contourf(rgrid, zgrid, log.(norm_factor * H_counter)', [-11,-6,-5.5,-5.0,-4.5,-4.0,-3.0,-2.0,-0.0], cmap="jet")
gca().invert_yaxis()

colorbar(cs,location="bottom", ax=ax, label = "log B")

display(gcf())
close(gcf())
## horizontal cross section
# Hz_cross = Hz[:,31]
# Hr_cross = Hr[:,31]
# thetagrid = (0:1:360)./360 * 2 * π

# Hfield_params = reduce(hcat, SNMRForward.co_counter_field.(Hz_cross, Hr_cross, 13*π/36, θ) for θ in thetagrid)
# Hco = first.(Hfield_params)
# ζ = last.(Hfield_params)
# H_counter = reshape([a[2] for a in Hfield_params[:]], size(Hco)...)

# fig = figure()
# ax = fig.add_subplot(projection="polar")
# μ_G = 4*π * 1e-3
# contourf(thetagrid .+ π/2, rgrid, real.(μ_G * Hco * 300) / 0.0319839, levels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1])
# display(gcf())

##
μ = SNMRForward.μ0
kernel = SNMRForward.point_kernel.(10, μ * Hco, μ * H_counter, ζ, ωl)

kernel *= SNMRForward.mag_factor(300.0) * ωl/SNMRForward.γh

# ## get the full perpendicular field
# Hperp = SNMRForward.perpendicular_field.(Hz,Hr,13*π/36,π/2)
# Hperp_z = imag.(first.(Hperp))

# figure()
# pcolor(rgrid, zgrid, Hperp_z', vmin = -0.06, vmax=0.06)
# gca().invert_yaxis()
# display(gcf())
# close(gcf())

# ##

# H_counter = reshape([a[2] for a in Hfield_params[:]], size(Hco)...)

# figure()
# pcolor(rgrid, zgrid, real.(H_counter)', vmin = 0, vmax=0.06)
# gca().invert_yaxis()
# display(gcf())
# close(gcf())

# ##
figure()
cs = contourf(rgrid, zgrid, real.(kernel)'/3.89619e-10, levels = [-0.1,-0.05,-0.02,-0.01,0.01,0.02,0.05,0.1],cmap = "gist_stern")
gca().invert_yaxis()
colorbar(cs)
display(gcf())

##
n_theta_points = 100
thetagrid = range(0, 2*pi, length=n_theta_points)

#radial integral scale
dr = (rgrid[2] - rgrid[1]) * rgrid
#azimuthal integral scale
dtheta = thetagrid[2] - thetagrid[1]

k1d = zeros(ComplexF64, size(zgrid)...)
# This should go in a "1D kernel" function later

q = 10
# ωl = 2*π*2.5e3
for (i_th, θ) = enumerate(thetagrid)
    Hparams = SNMRForward.co_counter_field.(Hz, Hr, ϕ, θ)

    Hco = first.(Hparams)
    ζ = last.(Hparams)
    Hctr = reshape([a[2] for a in Hparams[:]], size(Hco)...)

    kernel = SNMRForward.point_kernel.(q, μ * Hco, μ * Hctr, ζ, ωl)
    k1d += dtheta*transpose(kernel)*dr
end

##

Be = ωl/SNMRForward.γh
m0 = SNMRForward.mag_factor(300) * Be


##
full_kernel = k1d * m0 # in V/m
figure()
plot(real.(full_kernel * 10^9), zgrid)
gca().set_yscale("log")
gca().invert_yaxis()
display(gcf())

## tipping angle for comparison with literature
θ = -π/2
Hparams = SNMRForward.co_counter_field.(Hz, Hr, ϕ, θ)

Hco = first.(Hparams)
ζ = last.(Hparams)
Hctr = reshape([a[2] for a in Hparams[:]], size(Hco)...)

q = 10
αt = ((SNMRForward.γh * q * μ * Hco) .% (2*π))/π * 180

figure()
cs = contourf(rgrid, zgrid, αt', levels=[0,45,90,135,210,225,270,315,360])
gca().invert_yaxis()
xlabel("distance from loop centre (m)")
ylabel("depth (m)")
colorbar(cs, label = "tipping angle (degrees)")
display(gcf())
savefig("tipping.png")

## wrap it up in a function
function kernel_1d(q, ϕ, ωl, Hz, Hr)
    n_theta_points = 100
    thetagrid = range(0, 2*pi, length=n_theta_points)

    #radial integral scale
    dr = (rgrid[2] - rgrid[1]) * rgrid
    #azimuthal integral scale
    dtheta = thetagrid[2] - thetagrid[1]

    k1d = zeros(ComplexF64, size(Hz,2))

    for (i_th, θ) = enumerate(thetagrid)
        Hparams = SNMRForward.co_counter_field.(Hz, Hr, ϕ, θ)

        Hco = first.(Hparams)
        ζ = last.(Hparams)
        Hctr = reshape([a[2] for a in Hparams[:]], size(Hco)...)

        kernel = SNMRForward.point_kernel.(q, μ * Hco, μ * Hctr, ζ, ωl)
        k1d += dtheta*transpose(kernel)*dr
    end
    k1d
end

## contour plot of 1D kernel (cf. fig. 5.6, Hertrich)
qgrid = [0.1,0.25,0.5,0.75,1,2,3,4,5,6,7,8,9,10,11] .* 2
ϕ = 12*π/36
# ωl = 2*π*2.5e3

kq = reduce(hcat, kernel_1d(q, ϕ, ωl, Hz, Hr) for q in qgrid)

##
fig, ax = subplots(1,1,figsize=(5,10))
contourf(qgrid, zgrid, 10^9 * real.(kq*m0), cmap="RdBu_r", levels=[-150,-100,-50,-25,0,25,50,100,150])
gca().invert_yaxis()
xlabel("q (A s)")
ylabel("Depth (m)")
colorbar(label = "Real part of 1D kernel (nV/m)")
display(gcf())
close(gcf())

## "log sensitivity" (Fig. 5.5, Hertrich)
fig, ax = subplots(1,1,figsize=(7,10))
contourf(qgrid, zgrid, log10.(10^9 * abs.(kq*m0)), levels=[-1,0.1,0.6,1.15,1.7,2.5], cmap="jet")
gca().invert_yaxis()
xlabel("Pulse moment (A s)")
ylabel("Depth (m)")
colorbar(label = "log(1D sensitivity (nV/m))")
display(gcf())
close(gcf())

## do an actual forward model
dz = zgrid[2] - zgrid[1]
fwd_kernel = kq * m0
## 10 - 20 m
w = zeros(length(zgrid))
w[(zgrid .>= 10) .& (zgrid .<= 20)] .= 1

response = transpose(fwd_kernel) * w * dz

fig, ax = subplots(1,3, figsize=(15,5))
sca(ax[1])
plot(qgrid, real.(response))
title("Saturated layer 10 - 20 m")
ylabel("Response voltage (V)")
xlabel("Pulse moment (A s)")
display(gcf())
## 30 - 45 m
w = zeros(length(zgrid))
w[(zgrid .>= 30) .& (zgrid .<= 45)] .= 1

response = transpose(fwd_kernel) * w * dz

sca(ax[2])
plot(qgrid, real.(response))
title("Saturated layer 30 - 45 m")
xlabel("Pulse moment (A s)")
display(gcf())
## 60 - 80 m
w = zeros(length(zgrid))
w[(zgrid .>= 60) .& (zgrid .<= 80)] .= 1

response = transpose(fwd_kernel) * w * dz

sca(ax[3])
plot(qgrid, real.(response))
title("Saturated layer 60 - 80 m")
xlabel("Pulse moment (A s)")
display(gcf())
savefig("Forwards.png")
##

## horizontal cross section of the kernel
q = 10
Hz_cross = Hz[:,6]
Hr_cross = Hr[:,6]
thetagrid = (0:1:360)./360 * 2 * π

Hfield_params = reduce(hcat, SNMRForward.co_counter_field.(Hz_cross, Hr_cross, 13*π/36, θ) for θ in thetagrid)
Hco = first.(Hfield_params)
ζ = last.(Hfield_params)
H_counter = reshape([a[2] for a in Hfield_params[:]], size(Hco)...)
kernel_cross = m0 * SNMRForward.point_kernel.(q, μ*Hco, μ*H_counter, ζ, ωl) * 10^9

fig = figure()
ax = fig.add_subplot(projection="polar")
ax.set_theta_zero_location("N")
ax.set_rticks([25,50,75,100])
ax.set_yticklabels(["25","50","75","r = 100 m"])
pcolormesh(thetagrid, rgrid, real.(kernel_cross)/0.059512, cmap="RdBu", vmin=-1.0, vmax = 1.0)
colorbar(label="Normalised real kernel")
display(gcf())
## 1d kernels at q = 1, 5, 10, 15 As
qset = [1,5,10,15]

kset = [kernel_1d(q, ϕ, ωl, Hz, Hr) for q in qset]

##
fig, ax = subplots(1,length(kset), figsize = (15,5))

for (i, (q,kern)) = enumerate(zip(qset,kset))
    sca(ax[i])
    if i == 1
        ylabel("depth (m)")
    end
    xlabel("real kernel (nV/m)")
    title("q = $q As")
    full_kernel = kern*m0 * 10^9
    plot(real.(full_kernel), zgrid)
    gca().invert_yaxis()
    gca().set_yscale("log")
end
display(gcf())

## Actually use the defined structs in the package to do forward modelling
condLEM = SNMRForward.ConductivityModel(σ, d)
qgrid = [0.1,0.25,0.5,0.75,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
# ϕ = 2*π/3
ϕ = 13*π/36
F = SNMRForward.MRSForward(R, zgrid, qgrid, ϕ, Be, condLEM)
##
w = zeros(length(zgrid))
w[(zgrid .>= 30) .& (zgrid .<= 45)] .= 1
data = SNMRForward.forward(F,w)

figure()
plot(qgrid,real.(data))
display(gcf())

##