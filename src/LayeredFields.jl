## response function (computes "B-response"/impedance and alpha vals for layered earth)
RealVector = Vector{T} where T <: Real
function responses(σ::RealVector, d::RealVector, κ::RealVector, ωl::Real)
    B = zeros(ComplexF64, length(κ), length(σ));
    α = [sqrt(k^2 - im * ωl * μ0 * σm)  for k in κ, σm in σ]
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
# function phiz(phi0, Bresponse, α, d, zgrid)
#     #schelkunoff potential at the top of each layer
#     phi_tops = zeros(ComplexF64, size(Bresponse)...)
#     phi_tops[:,1] = phi0
#     for m=1:size(phi_tops,2)-1
#         phi_tops[:,m+1] = phi_tops[:,m] .* (α[:,m] + Bresponse[:,m])./(α[:,m] + Bresponse[:,m+1]) .* exp.(-α[:,m]*d[m])
#     end
#     #convert depth to z-values for the bottom of each layer
#     z_interface = cumsum(d)

#     phiz = zeros(ComplexF64, length(phi0), length(zgrid))
#     phipz = zeros(ComplexF64, length(phi0), length(zgrid))

#     layer_m = 1
#     coeff_minus = phi_tops[:,1] .* (1 .+ Bresponse[:,1]./α[:,1])/2
#     coeff_plus = 0
#     if length(d) > 1
#         coeff_plus = phi_tops[:,2] .* (1 .- Bresponse[:,2]./α[:,2])/2
#     end
#     h_m = 0
#     h_mp1 = z_interface[1]

#     for (iz, z) = enumerate(zgrid)
#         if z > h_mp1
#             coeff_minus = phi_tops[:,layer_m] .* (1 .+ Bresponse[:,layer_m]./α[:,layer_m])/2
#             coeff_plus = 0
#             if length(d) > layer_m
#                 coeff_plus = phi_tops[:,layer_m+1] .* (1 .- Bresponse[:,layer_m+1]./α[:,layer_m+1])/2
#             end
#             h_m = h_mp1
#             h_mp1 = z_interface[layer_m]
#         end
        
#         phiz[:,iz] = coeff_minus .* exp.(-α[:,layer_m]*(z-h_m)) +
#                         coeff_plus .* exp.(α[:,layer_m]*(z-h_mp1))
#         phipz[:,iz] = α[:,layer_m].* (coeff_plus .* exp.(α[:,layer_m]*(z-h_mp1))
#                         - coeff_minus .* exp.(-α[:,layer_m]*(z-h_m)))
#     end

#     phiz, phipz
# end

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
    z_interface = [cumsum(d); Inf]

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

#spatial-domain magnetic H fields in vertical and radial directions for loop radius R
function magfields(R::Real, ω::Real, σ::AbstractVector{<:Real}, d::AbstractVector{<:Real}, rgrid::AbstractVector{<:Real}, zgrid::AbstractVector{<:Real})
    #define k grid for j0 and j1 kernels
    #j0        
    kj0 = reduce(hcat, Filter_base_801/r for r in rgrid)
    rj0 = rgrid' .* ones(length(Filter_base_801), length(rgrid))
    #j1
    kj1 = Filter_base_801 / R

    ## determine whether to use j0 or j1 for hz at each r
    use_j1 = [r < R for r in rgrid]
    j1_inds = cumsum(use_j1)

    ## compute layer responses to propagate fields
    Bj0, αj0 = responses(σ, d, kj0[:], ω)
    Bj1, αj1 = responses(σ, d, kj1, ω)

    ## propagation coefficients for phi and phi' (Hz and Hr)
    phicj0, phipj0 = phi_coeffs(Bj0, αj0, d, zgrid)
    phicj1, phipj1 = phi_coeffs(Bj1, αj1, d, zgrid)

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

    # Hz_kernel_j1 =  lowpass_k.(kj1,R) .* kj1 .^ 3 .* permutedims(phiz_j1, (1,3,2))
    Hz_kernel_j1 = kj1 .^ 3 .* permutedims(phiz_j1, (1,3,2))

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
    
    (Hz,Hr)
end