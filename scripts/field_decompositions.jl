using SNMRForward

R = 10
rgrid = 0.1:0.2:1.5*R
zgrid = 0:0.1:5*R

## conductive half-space, at 2 kHz
ωl = 2.0e3 #Hz, typical for Earth's field strength
d = [Inf]
σ = [5]

#define k grid for j0 and j1 kernels
k_vals = reduce(hcat, Filter_base / r for r in rgrid)
k_vals_j1 = Filter_base / R

#k space potentials and derivatives for each r value


##
B_halfspace, alpha_halfspace = responses(σ, d, k_vals, ωl)