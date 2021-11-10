using Revise, SNMRForward

R = 50
rgrid = 0.1:0.5:2*R
zgrid = 0.1*R:0.5:2*R

## conductive half-space, at 2 kHz
ωl = 2*π*2.5e3 #Hz, typical for Earth's field strength
d = [Inf]
σ = [0.05]

#define k grid for j0 and j1 kernels

#j0        
kj0 = reduce(hcat, SNMRForward.Filter_base_801/r for r in rgrid)
rj0 = rgrid' .* ones(length(SNMRForward.Filter_base_801), length(rgrid))
#j1
kj1 = SNMRForward.Filter_base_801 / R

## determine whether to use j0 or j1 for hz at each r
use_j1 = [r < R for r in rgrid]
j1_inds = cumsum(use_j1)

## compute layer responses to propagate fields
Bj0, αj0 = SNMRForward.responses(σ, d, kj0[:], ωl)
Bj1, αj1 = SNMRForward.responses(σ, d, kj1, ωl)

## propagation coefficients for phi and phi' (Hz and Hr)
phicj0, phipj0 = SNMRForward.phi_coeffs(Bj0, αj0, d, zgrid)
phicj1, phipj1 = SNMRForward.phi_coeffs(Bj1, αj1, d, zgrid)

## calculate z = 0 for j0 and j1 kernels

phif_j0 = SNMRForward.phi_free.(kj0[:], 0, R)

phi0_j0 = 2 * kj0[:]./(kj0[:] .+ Bj0) .* phif_j0

##
rj1 = rgrid[use_j1]
##
phif_j1 = reduce(hcat, SNMRForward.phi_free_j1k.(kj1, 0, R, r) for r in rj1)

phi0_j1 = 2 * kj1 ./ (kj1 .+ Bj1) .* phif_j1

## propagate through halfspace
phiz_j0 = phi0_j0 .* phicj0
phipz_j0 = phi0_j0 .* phipj0

##
phiz_j1 = cat([phi0_j1[:,ir] .* phicj1 for ir=1:size(phi0_j1,2)]...; dims=3)

##low-pass filters
filterspacing = SNMRForward.Filter_base_801[2] / SNMRForward.Filter_base_801[1]
kcut(r) = pi/((filterspacing-1)*r)
lowpass_k(k,r) = 1/(1 + (k/kcut(r))^3)

Hz_kernel_j1 =  lowpass_k.(kj1,R) .* kj1 .^ 3 .* permutedims(phiz_j1, (1,3,2))


## compute Hz and Hr field kernels for each k value
Hz_kernel = reshape(kj0[:].^3 .* phiz_j0, size(kj0)..., length(zgrid))
# note the Hr kernel is actually in j1 space
# due to the hankel transform property used to derive it
Hr_kernel = reshape(-kj0[:].^2 .* phipz_j0, size(kj0)..., length(zgrid))

# recover z and r fields on z and r grid
Hz = zeros(ComplexF64, length(rgrid), length(zgrid))
Hr = zeros(ComplexF64, length(rgrid), length(zgrid))
for ir in 1:length(rgrid), iz in 1:length(zgrid)
    if use_j1[ir]
        ir1 = j1_inds[ir]
        Hz[ir,iz] = 1/R * SNMRForward.Filter_J1_801' * Hz_kernel_j1[:,ir1,iz]
    else
        Hz[ir, iz] = 1/rgrid[ir] * SNMRForward.Filter_J0_801' * Hz_kernel[:,ir,iz]
    end
    Hr[ir, iz] = 1/rgrid[ir] * SNMRForward.Filter_J1_801' * Hr_kernel[:,ir,iz]
end

## plot the real part of the fields
using PyPlot
fig, ax = subplots(1,2)
sca(ax[1])
pcolor(rgrid, zgrid, real.(Hz)', vmin = -0.06, vmax=0.06)
gca().invert_yaxis()
sca(ax[2])
pcolor(rgrid, zgrid, real.(Hr)', vmin = -0.06, vmax=0.06)
gca().invert_yaxis()
display(gcf())
close(gcf())
## say the mag. field inclination is 60 degrees. calculate corotating part and phase lag
Hfield_params = SNMRForward.co_counter_field.(Hz, Hr, 13*π/36, -π/2)

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

##
maximum(imag.(Hco))

##
μ = SNMRForward.mu_0
kernel = SNMRForward.point_kernel.(10, μ * Hco, μ * H_counter, ζ, ωl)

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
pcolor(rgrid, zgrid, imag.(kernel)')
gca().invert_yaxis()
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

ϕ = 13*π/36
q = 3
ωl = 2.5e3
for (i_th, θ) = enumerate(thetagrid)
    Hparams = SNMRForward.co_counter_field.(Hz, Hr, ϕ, θ)

    Hco = first.(Hparams)
    ζ = last.(Hparams)
    Hctr = reshape([a[2] for a in Hfield_params[:]], size(Hco)...)

    kernel = SNMRForward.point_kernel.(q, μ * Hco, μ * Hctr, ζ, ωl)
    k1d += dtheta*kernel'*dr
end
figure()
plot(zgrid,abs.(k1d))
display(gcf())