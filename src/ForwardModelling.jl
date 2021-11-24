#1D forward modelling object for
#MRS experiment with circular loop
struct MRSForward
    kernel::Matrix{<:Number}      #kernel (precomputed and stored)
    qgrid::AbstractVector{<:Real}
    zgrid::AbstractVector{<:Real}         #loop radius
end

struct ConductivityModel
    σ::Vector{<:Real}
    d::Vector{<:Real}
    ConductivityModel(σ,t) = (length(t) != length(σ) - 1 ?
            error("thickness vector must be one shorter than conductivity vector") :
            new(σ,t))
end

function MRSForward(R::Real, zgrid::AbstractVector{<:Real},
    qgrid::AbstractVector{<:Real}, ϕ::Real, Be::Real,
    condLEM::ConductivityModel; nrvals = 200, temp=300.0)
    #radial evaluation values
    R <= 0 && error("loop radius must be positive")
    rmin = 0.01*R
    rmax = 4*R
    rgrid = range(rmin, stop=rmax, length=nrvals)
    m0 = mag_factor(temp) * Be
    ωl = γh * Be
    (Hz, Hr) = magfields(R, ωl, condLEM.σ, condLEM.d, rgrid, zgrid)
    kq = reduce(hcat, kernel_1d(q, ϕ, ωl, Hz, Hr, rgrid) for q in qgrid)
    dz = zgrid[2:end] - zgrid[1:end-1]
    k_mat = m0 * 0.5 * dz .* (kq[1:end-1, :] .+ kq[2:end,:])
    MRSForward(k_mat, qgrid, zgrid)
end

function forward(F::MRSForward, m::Vector{<:Real})
    #trapezoid integration of the kernel * the model vector
    #volume element is included in the stored kernel
    m_trap = 0.5*(m[1:end-1] .+ m[2:end])
    transpose(F.kernel)*m_trap
end

