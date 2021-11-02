using SNMRForward

R = 10
rgrid = 0.1:0.2:1.5*R
zgrid = 0:0.1:1.5*R

## conductive half-space, at 2 kHz
ωl = 2.0e3 #Hz, typical for Earth's field strength
d = [Inf]
σ = [5]

#define k grid for j0 and j1 kernels
k_vals = reduce(hcat, SNMRForward.Filter_base / r for r in rgrid)
k_vals_j1 = SNMRForward.Filter_base / R



## compute layer responses to propagate fields
B_halfspace, α_halfspace = SNMRForward.responses(σ, d, k_vals[:], ωl)
B_halfspace_j1, α_halfspace_j1 = SNMRForward.responses(σ, d, k_vals_j1, ωl)

## 
phif_j0 = SNMRForward.phi_free.(k_vals[:], 0, R)

phi0_j0 = 2 * k_vals[:]./(k_vals[:] .+ B_halfspace[:,1]) .* phif_j0


## propagate through halfspace
phiz_j0, phipz_j0 = SNMRForward.phiz(phi0_j0, B_halfspace, α_halfspace, d, zgrid)

## compute Hz and Hr field kernels for each k value
Hz_kernel = reshape(k_vals[:].^3 .* phiz_j0, size(k_vals)..., length(zgrid))
# note the Hr kernel is actually in j1 space
# due to the hankel transform property used to derive it
Hr_kernel = reshape(- k_vals[:].^2 .* phipz_j0, size(k_vals)..., length(zgrid))

# recover z and r fields on z and r grid
Hz = zeros(ComplexF64, length(rgrid), length(zgrid))
Hr = zeros(ComplexF64, length(rgrid), length(zgrid))
for ir in 1:length(rgrid), iz in 1:length(zgrid)
    Hz[ir, iz] = 1/rgrid[ir] * SNMRForward.Filter_J0' * Hz_kernel[:,ir,iz]
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
Hfield_params = SNMRForward.co_counter_field.(Hz, Hr, π/3, 0)

Hco = first.(Hfield_params)
ζ = last.(Hfield_params)

##
fig, ax = subplots(1,2)
sca(ax[1])
pcolor(rgrid, zgrid, real.(Hco)', vmin = -0.06, vmax=0.06)
gca().invert_yaxis()
display(gcf())
sca(ax[2])
pcolor(rgrid, zgrid, ζ', vmin = -π, vmax=π)
gca().invert_yaxis()
display(gcf())
close(gcf())

