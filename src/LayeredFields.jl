## response function (computes "B-response"/impedance and alpha vals for layered earth)
RealVector = Vector{T} where T <: Real
function responses(σ::RealVector, d::RealVector, κ::RealVector, ωl::Real)
    B = zeros(ComplexF64, length(κ), length(σ));
    α = [sqrt(k^2 - im * ωl * mu_0 * σm)  for k in κ, σm in σ]
    B[:,end] = α[:,end]
    m = size(B,2) - 1

    while m > 0
        thamdm = tanh.(α[:,m]*d[m])
        B[:,m] = α[:,m] .* (B[:,m+1] .+ α[:,m].*thamdm)./(α[:,m] .+ B[:,m+1].*thamdm)
        m -= 1
    end
    return B, α
end


#Schelkunoff potential and its derivative at each z
function phiz(phi0, Bresponse, α, d, zgrid)
    #schelkunoff potential at the top of each layer
    phi_tops = zeros(ComplexF64, size(Bresponse)...)
    phi_tops[:,1] = phi0
    for m=1:size(phi_tops,2)-1
        phi_tops[:,m+1] = phi_tops[:,m] .* (α[:,m] + Bresponse[:,m])./(α[:,m] + Bresponse[:,m+1]) .* exp.(-α[:,m]*d[m])
    end
    #convert depth to z-values for the bottom of each layer
    z_interface = cumsum(d)

    phiz = zeros(ComplexF64, length(phi0), length(zgrid))
    phipz = zeros(ComplexF64, length(phi0), length(zgrid))

    layer_m = 1
    coeff_minus = phi_tops[:,1] .* (1 .+ Bresponse[:,1]./α[:,1])/2
    coeff_plus = 0
    if length(d) > 1
        coeff_plus = phi_tops[:,2] .* (1 .- Bresponse[:,2]./α[:,2])/2
    end
    h_m = 0
    h_mp1 = z_interface[1]

    for (iz, z) = enumerate(zgrid)
        if z > h_mp1
            coeff_minus = phi_tops[:,layer_m] .* (1 .+ Bresponse[:,layer_m]./α[:,layer_m])/2
            coeff_plus = 0
            if length(d) > layer_m
                coeff_plus = phi_tops[:,layer_m+1] .* (1 .- Bresponse[:,layer_m+1]./α[:,layer_m+1])/2
            end
            h_m = h_mp1
            h_mp1 = z_interface[layer_m]
        end
        
        phiz[:,iz] = coeff_minus .* exp.(-α[:,layer_m]*(z-h_m)) +
                        coeff_plus .* exp.(α[:,layer_m]*(z-h_mp1))
        phipz[:,iz] = α[:,layer_m].* (coeff_plus .* exp.(α[:,layer_m]*(z-h_mp1))
                        - coeff_minus .* exp.(-α[:,layer_m]*(z-h_m)))
    end

    phiz, phipz
end

#propagation coefficients in a layered earth
#multiply these by field values at z = 0
function phi_coeffs(Bresponse, α, d, zgrid)
    #schelkunoff potential at the top of each layer
    phi_tops = zeros(ComplexF64, size(Bresponse)...)
    phi_tops[:,1] .= 1
    for m=1:size(phi_tops,2)-1
        phi_tops[:,m+1] = phi_tops[:,m] .* (α[:,m] + Bresponse[:,m])./(α[:,m] + Bresponse[:,m+1]) .* exp.(-α[:,m]*d[m])
    end
    #convert depth to z-values for the bottom of each layer
    z_interface = cumsum(d)

    phiz = zeros(ComplexF64, size(Bresponse,1), length(zgrid))
    phipz = zeros(ComplexF64, size(Bresponse,1), length(zgrid))

    layer_m = 1
    coeff_minus = phi_tops[:,1] .* (1 .+ Bresponse[:,1]./α[:,1])/2
    coeff_plus = 0
    if length(d) > 1
        coeff_plus = phi_tops[:,2] .* (1 .- Bresponse[:,2]./α[:,2])/2
    end
    h_m = 0
    h_mp1 = z_interface[1]

    for (iz, z) = enumerate(zgrid)
        if z > h_mp1
            coeff_minus = phi_tops[:,layer_m] .* (1 .+ Bresponse[:,layer_m]./α[:,layer_m])/2
            coeff_plus = 0
            if length(d) > layer_m
                coeff_plus = phi_tops[:,layer_m+1] .* (1 .- Bresponse[:,layer_m+1]./α[:,layer_m+1])/2
            end
            h_m = h_mp1
            h_mp1 = z_interface[layer_m]
        end
        
        phiz[:,iz] = coeff_minus .* exp.(-α[:,layer_m]*(z-h_m)) +
                        coeff_plus .* exp.(α[:,layer_m]*(z-h_mp1))
        phipz[:,iz] = α[:,layer_m].* (coeff_plus .* exp.(α[:,layer_m]*(z-h_mp1))
                        - coeff_minus .* exp.(-α[:,layer_m]*(z-h_m)))
    end

    phiz, phipz
end