#kernel function evaluated at a point, given co- and counter-rotating fields
#and phase lag
#must be scaled by equilibrium magnetisation (depends on magnitude of earth field)
point_kernel(q, Bplus, Bminus, ζ, ωl) = 
    -2 * ωl * sin(γh * q * Bplus) * Bminus * exp(-2*im*ζ)