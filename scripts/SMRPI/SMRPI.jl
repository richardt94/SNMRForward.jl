# wrapper to integrate SMRPI 1D forward model with
# HiQGA/transD_GP
module SMRPI
using transD_GP.AbstractOperator, transD_GP.CommonToAll
import transD_GP.AbstractOperator.get_misfit
import transD_GP.Model, transD_GP.Options
import transD_GP.ModelNuisance, transD_GP.OptionsNuisance
using SNMRForward, Random

export γh, ConductivityModel, MRSForward, MRSForward_square
export newSMRSounding, create_synthetic

abstract type SMRSounding <: Operator1D end

mutable struct SMRSoundingKnown <: SMRSounding
    V0 :: Vector{<:Real} #sounding curve
    ϕ :: Vector{<:Real} #phases
    σ_V0 :: Vector{<:Real}
    σ_ϕ :: Vector{<:Real}
    Fm :: SNMRForward.MRSForward
    linearsat :: Bool
    amponly :: Bool
end

mutable struct SMRSoundingUnknown <: SMRSounding
    # this struct doesn't contain noise variances
    # - a maximum likelihood estimate is computed by the sampler instead
    V0 :: Vector{<:Real}
    ϕ :: Vector{<:Real}
    Fm :: SNMRForward.MRSForward
    σd_diag :: Vector{<:Real}
    linearsat :: Bool
    amponly :: Bool
end

function newSMRSounding(V0, ϕ, Fm; σ_V0=nothing, σ_ϕ=nothing, mult=false, linearsat=false, amponly=false)
    (length(V0) != length(ϕ) && 
    throw(ArgumentError("V0 and ϕ must have same length")))
    if !isnothing(σ_V0) || !isnothing(σ_ϕ)
        if isnothing(σ_V0) || isnothing(σ_ϕ)
            throw(ArgumentError("σ_V0 and σ_ϕ must both be provided, or neither."))
        end
        if length(V0) != length(σ_V0) || length(ϕ) != length(σ_ϕ)
            throw(ArgumentError("σ_V0 and σ_ϕ must have the same length as the associated data"))
        end
        return SMRSoundingKnown(V0, ϕ, σ_V0, σ_ϕ, Fm, linearsat, amponly)
    end

    if mult
        σd_diag = [V0; ones(length(ϕ))]
    else
        σd_diag = [ones(length(V0)); 1 ./ V0]
    end
    SMRSoundingUnknown(V0, ϕ, Fm, σd_diag, linearsat, amponly)
end

function get_misfit(m::Model, opt::Options, S::SMRSounding)
    opt.debug && return 0.0
    get_misfit(S.linearsat ? m.fstar[:] : 10 .^ m.fstar[:], S)
end
# above defined function and type signature MUST be defined

function get_misfit(m::Model, mn::ModelNuisance, opt::Union{Options,OptionsNuisance}, S::SMRSounding)
    #the "nuisance" in this case is a constant offset phase
    opt.debug && return 0.0
    offset_ϕ = mn.nuisance[1]
    get_misfit(S.linearsat ? m.fstar[:] : 10 .^ m.fstar[:], S, offset_ϕ = offset_ϕ)
end

function get_misfit(w::Vector{<:Real}, S::SMRSounding; offset_ϕ = 0.)
    response = SNMRForward.forward(S.Fm, w)
    Vres = abs.(response)
    ϕres = angle.(response) .+ offset_ϕ
    if isa(S, SMRSoundingKnown)
        #use provided noise
        residual = [(S.V0 .- Vres)./S.σ_V0]
        S.amponly || (residual = [residual; rem2pi.(S.ϕ - ϕres, RoundNearest)./S.σ_ϕ])
        return residual' * residual / 2
    else
        #maximum-likelihood estimate of multiplicative noise
        residual = [S.V0 - Vres]./S.σd_diag[1:length(Vres)]
        S.amponly || (residual = [residual; rem2pi.(S.ϕ - ϕres, RoundNearest)./S.σd_diag[length(Vres)+1:end]])
        return length(residual)/2 * log(residual' * residual)
    end
end

function create_synthetic(w::Vector{<:Real}, σ::Vector{<:Real}, t::Vector{<:Real},
            Be::Real, ϕ::Real, R::Real, zgrid::Vector{<:Real}, qgrid::Vector{<:Real}
    ; noise_frac = 0.05, θ = 0., square=false, noise_mle = false, mult = false, linearsat=false,
    amponly=false, offset_ϕ = 0.)
    ct = SNMRForward.ConductivityModel(σ, t)

    F = (square ?
        SNMRForward.MRSForward(R, zgrid, qgrid, ϕ, Be, ct) :
        SNMRForward.MRSForward_square(R, zgrid, qgrid, ϕ, θ, Be, ct))

    synth_data = SNMRForward.forward(F,w)
    if mult
        σ_V0 = noise_frac * abs.(synth_data)
        σ_ϕ = noise_frac * ones(size(synth_data))
    else
        σ = noise_frac * maximum(abs.(synth_data))
        σ_V0 = σ * ones(size(synth_data))
        σ_ϕ = σ ./ abs.(synth_data)
    end
    noisy_V0 = abs.(synth_data) .+ σ_V0 .* randn(size(abs.(synth_data)))
    noisy_ϕ = angle.(synth_data) .+ σ_ϕ .* randn(size(abs.(synth_data))) .+ offset_ϕ

    if noise_mle
        σ_V0 = nothing
        σ_ϕ = nothing
    end

    newSMRSounding(noisy_V0, noisy_ϕ, F, σ_V0=σ_V0, σ_ϕ=σ_ϕ, mult=mult, linearsat=linearsat, amponly=amponly)
end


end