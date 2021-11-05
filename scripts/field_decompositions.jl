using Revise, SNMRForward

R = 50
rgrid = 0.1:0.5:2*R
zgrid = 0:0.5:2*R

## conductive half-space, at 2 kHz
ωl = 2.5e3 #Hz, typical for Earth's field strength
d = [Inf]
σ = [0.05]

#define k grid for j0 and j1 kernels

#j0        
kj0 = reduce(hcat, SNMRForward.Filter_base/r for r in rgrid)
rj0 = rgrid' .* ones(length(SNMRForward.Filter_base), length(rgrid))
#j1
kj1 = SNMRForward.Filter_base / R

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
filterspacing = SNMRForward.Filter_base[2] / SNMRForward.Filter_base[1]
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
        Hz[ir,iz] = 1/R * SNMRForward.Filter_J1' * Hz_kernel_j1[:,ir1,iz]
    else
        Hz[ir, iz] = 1/rgrid[ir] * SNMRForward.Filter_J0' * Hz_kernel[:,ir,iz]
    end
    Hr[ir, iz] = 1/rgrid[ir] * SNMRForward.Filter_J1' * Hr_kernel[:,ir,iz]
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
Hfield_params = SNMRForward.co_counter_field.(Hz, Hr, 13*π/36, π/2)

Hco = first.(Hfield_params)
ζ = real.(last.(Hfield_params))


H_counter = reshape([a[2] for a in Hfield_params[:]], size(Hco)...)

##
fig, ax = subplots(1,2)
sca(ax[1])
pcolor(rgrid, zgrid, log.(real.(Hco))', vmin = -7.5, vmax=-2)
gca().invert_yaxis()
display(gcf())
sca(ax[2])
pcolor(rgrid, zgrid, log.(real.(H_counter))', vmin = -7.5, vmax=-2)
gca().invert_yaxis()
display(gcf())
close(gcf())

##

figure()
pcolor(rgrid, zgrid, ζ', vmin = -π/2, vmax=π/2)
gca().invert_yaxis()
display(gcf())
##
maximum(imag.(Hco))



## get the full perpendicular field
Hperp = SNMRForward.perpendicular_field.(Hz,Hr,13*π/36,π/2)
Hperp_z = imag.(first.(Hperp))

figure()
pcolor(rgrid, zgrid, Hperp_z', vmin = -0.06, vmax=0.06)
gca().invert_yaxis()
display(gcf())
close(gcf())

##

H_counter = reshape([a[2] for a in Hfield_params[:]], size(Hco)...)

figure()
pcolor(rgrid, zgrid, real.(H_counter)', vmin = 0, vmax=0.06)
gca().invert_yaxis()
display(gcf())
close(gcf())

##