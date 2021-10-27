using SpecialFunctions

#schelkunoff potential from a unit current loop in free space
#with radius R
#defined in the frequency-horizontal wavenumber domain
#for z > 0
function phi_free(κ, z, R)
    (R * besselj1(κ * R) ./ (2 * κ^2)) * exp(-κ*z)
end