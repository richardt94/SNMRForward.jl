# find components of mag. field perpendicular
# to earth's field, given inclination φ
# and azimuthal coordinate θ
# for cylindrically-symmetric fields with no azimuthal component
# (e.g. circular loop source)
function perpendicular_field(Bz, Br, φ, θ)
    #y is defined as magnetic west.
    #x is defined as
    #cross product of y and the unit earth field vector
    #(perpendicular to the earth field in the north-down plane)
    Bperp_x = -Bz * cos(φ) + Br * cos(θ) * sin(φ)
    Bperp_y = Br * sin(θ) 
    [Bperp_x; Bperp_y]
end

#for fields expressed in cartesian coordinates
#(e.g. square loop source)
#θ is angle between magnetic north and x axis
#ϕ is inclination
function perpendicular_field(Bx, By, Bz, ϕ, θ)
    Bperp_x = Bx * cos(θ) * sin(ϕ) + By * sin(θ) * sin(ϕ) + Bz * cos(ϕ)
    Bperp_y = By * cos(θ) - Bx * sin(θ)
    [Bperp_x; Bperp_y]
end

# find co-rotating and counter-rotating field parameters
# need |B+|, |B-| and the phase lag ζt
function co_counter_field(Bperp)
    Bstar = conj(Bperp)
    BdotB = transpose(Bperp) * Bperp
    BdotBstar =  transpose(Bperp) * Bstar
    BcrossBstar =  imag(Bperp[1] * Bstar[2] - Bperp[2] * Bstar[1])

    ζ = -real(angle(BdotB/abs(BdotB))/2)

    αt = real(sqrt((BdotBstar + abs(BdotB))/2))
    βt = real(sqrt((BdotBstar - abs(BdotB))/2))
    
    if BcrossBstar < 0
        βt *= -1
    end
    
    Bplus = (αt - βt)/2
    Bminus = (αt + βt)/2

    [Bplus; Bminus; ζ]
end

#convenience wrappers around the corotating/counter-rotating decomp
#for both cylindrical and cartesian coordinates fields
function co_counter_field(Bz, Br, φ, θ)
    Bperp = perpendicular_field(Bz, Br, φ, θ)
    co_counter_field(Bperp)
end

function co_counter_field(Bx, By, Bz, ϕ, θ)
    Bperp = perpendicular_field(Bx, By, Bz, ϕ, θ)
    co_counter_field(Bperp) 
end
