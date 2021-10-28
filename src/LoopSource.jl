using SpecialFunctions

#schelkunoff potential from a unit current loop in free space
#with radius R
#defined in the frequency-horizontal wavenumber domain
#for z > 0
function phi_free(κ, z, R)
    (R * besselj1(κ * R) ./ (2 * κ^2)) * exp(-κ*z)
end

#for inverse hankel transform using a J1 digital filter
#instead of a J0
function phi_free_j1k(κ, z, R, r)
    (R * besselj0(κ * r) ./ (2 * κ^2)) * exp(-κ*z)
end