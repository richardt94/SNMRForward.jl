# wrapper to integrate SMRPI 1D forward model with
# HiQGA/transD_GP
module SMRPI
using transD_GP.AbstractOperator, transD_GP.CommonToAll
import transD_GP.AbstractOperator.get_misfit
import transD_GP.Model, transD_GP.Options
import transD_GP.ModelNuisance, transD_GP.OptionsNuisance
using SNMRForward, Random, PyPlot

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

function newSMRSounding(V0, ϕ, Fm; σ_V0=nothing, σ_ϕ=nothing, mult=false, linearsat=false, amponly=false, showplot=false)
    (length(V0) != length(ϕ) && 
    throw(ArgumentError("V0 and ϕ must have same length")))
    if !isnothing(σ_V0) || !isnothing(σ_ϕ)
        if isnothing(σ_V0) || isnothing(σ_ϕ)
            throw(ArgumentError("σ_V0 and σ_ϕ must both be provided, or neither."))
        end
        if length(V0) != length(σ_V0) || length(ϕ) != length(σ_ϕ)
            throw(ArgumentError("σ_V0 and σ_ϕ must have the same length as the associated data"))
        end
        S = SMRSoundingKnown(V0, ϕ, σ_V0, σ_ϕ, Fm, linearsat, amponly)
        return S
    end

    if mult
        σd_diag = [V0; ones(length(ϕ))]
    else
        σd_diag = [ones(length(V0)); 1 ./ V0]
    end
    S = SMRSoundingUnknown(V0, ϕ, Fm, σd_diag, linearsat, amponly)
    showplot && plotdata(S)
    S
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
        residual = (S.V0 .- Vres)./S.σ_V0
        S.amponly || (residual = [residual; rem2pi.(S.ϕ - ϕres, RoundNearest)./S.σ_ϕ])
        return residual' * residual / 2
    else
        #maximum-likelihood estimate of multiplicative noise
        residual = (S.V0 - Vres)./S.σd_diag[1:length(Vres)]
        S.amponly || (residual = [residual; rem2pi.(S.ϕ - ϕres, RoundNearest)./S.σd_diag[length(Vres)+1:end]])
        return length(residual)/2 * log(residual' * residual)
    end
end

function create_synthetic(w::Vector{<:Real}, σ::Vector{<:Real}, t::Vector{<:Real},
            Be::Real, ϕ::Real, R::Real, zgrid::Vector{<:Real}, qgrid::Vector{<:Real}
    ; noise_frac = 0.05, θ = 0., square=false, noise_mle = false, mult = false, linearsat=false,
    amponly=false, offset_ϕ = 0., showplot=true, rseed=131)
    Random.seed!(rseed)
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
    S = newSMRSounding(noisy_V0, noisy_ϕ, F, σ_V0=σ_V0, σ_ϕ=σ_ϕ, mult=mult, linearsat=linearsat, amponly=amponly)
    if showplot 
        fig = plotmodelcurve(ct.σ, ct.d, w, zgrid, abs.(synth_data), angle.(synth_data), qgrid)
        plotdata(S, fig, iaxis=3, writelabel=false, msize=8)
    end
    S
end

## plotting stuff

function plotdata(S::SMRSounding; gridalpha=0.5, figsize=(6,3), msize=8)
    fig, ax = plt.subplots(1, 2, sharex=true, figsize=figsize)
    plotdata(S, fig, gridalpha=gridalpha, msize=msize)
    fig 
end

function plotdata(S::SMRSounding, fig; iaxis=1, gridalpha=0.5, writelabel=true, msize=8)
    ax = fig.axes
    ϕ = rem2pi.(S.ϕ, RoundNearest)
    if isa(S, SMRSoundingKnown)
        σ_ϕ = rem2pi.(S.σ_ϕ, RoundNearest)
        ax[iaxis].errorbar(S.Fm.qgrid, S.V0, S.σ_V0, linestyle="none", marker=".", elinewidth=1, capsize=3)
        ax[iaxis+1].errorbar(S.Fm.qgrid, ϕ, σ_ϕ, linestyle="none", marker=".", elinewidth=1, capsize=3)
    else
        ax[iaxis].plot(S.Fm.qgrid, S.V0, linestyle="none", marker="o", markersize=msize)
        ax[iaxis].plot(S.Fm.qgrid, S.V0, linestyle="none", marker=".", markersize=msize/2)
        ax[iaxis+1].plot(S.Fm.qgrid, ϕ, linestyle="none", marker="o", markersize=msize)
        ax[iaxis+1].plot(S.Fm.qgrid, ϕ, linestyle="none", marker=".", markersize=msize/2)
    end
    if writelabel # useful on the last iteration of model plotting
        writelabels(ax, iaxis, gridalpha)
        ax[iaxis+1].set_ylim(-pi, pi)
    end 
    fig.tight_layout()
end    

function plotmodelcurve(c, t, w, z, V0, ϕ, q; gridalpha=0.5, modelalpha=0.5,figsize=(10,3),
    lcolor="nocolor", writelabel=true)
    # conductivity and saturation with depth initialize
    fig = figure(figsize=(figsize))
    s1 = subplot(141)
    s2 = subplot(142, sharey=s1)
    s3 = subplot(143)
    s4 = subplot(144, sharex=s3)
    plotmodelcurve(c, t, w, z, V0, ϕ, q, fig, gridalpha=gridalpha, modelalpha=modelalpha, writelabel=writelabel, lcolor=lcolor)
    fig
end

function plotmodelcurve(c, t, w, z, V0, ϕ, q, fig; gridalpha=0.5, modelalpha=0.5, writelabel=true, lcolor="nocolor")
    # conductivity and saturation into axis
    ax = fig.axes
    zfromt = [0.; cumsum(t)]
    isempty(t) && push!(zfromt, maximum(z))
    ax[1].step([c;c[end]], zfromt)
    if writelabel
        ax[1].grid(b=true, which="both", alpha=gridalpha)
        ax[1].set_xlabel("conductivity S/m")
        ax[1].set_ylim(reverse(extrema(z)))
    end
    plotmodelcurve(w, z, V0, ϕ, q, fig; iaxis=2, gridalpha=gridalpha, modelalpha=modelalpha, writelabel=writelabel, lcolor=lcolor,)
end    

function plotmodelcurve(w, z, V0, ϕ, q; gridalpha=0.5, modelalpha=0.5,figsize=(10,3),
    lcolor="nocolor", writelabel=true)
    # saturation with depth initialize
    fig = figure(figsize=(figsize))
    s1 = subplot(131)
    s2 = subplot(132)
    s3 = subplot(133, sharex=s2)
    plotmodelcurve(w, z, V0, ϕ, q, fig, gridalpha=gridalpha, modelalpha=modelalpha, writelabel=writelabel, lcolor=lcolor)
    fig
end

function plotmodelcurve(w, z, V0, ϕ, q, fig; iaxis=1, gridalpha=0.5, modelalpha=0.5, writelabel=true, lcolor="nocolor")
    # saturation with depth into axis
    ax = fig.axes
    ax[iaxis].step(w, z, color=lcolor, alpha=modelalpha)
    if writelabel
        ax[iaxis].grid(b=true, which="both", alpha=gridalpha)
        ax[iaxis].set_xlabel("saturation")
        ax[iaxis].set_ylabel("depth m")
        ax[iaxis].set_ylim(reverse(extrema(z)))
    end    
    plotcurve(V0, ϕ, q, fig, iaxis=iaxis+1, gridalpha=gridalpha, modelalpha=modelalpha,
                    lcolor=lcolor, writelabel=writelabel)
end

function plotcurve(V0, ϕ, q, gridalpha=0.5)
    # sometimes you only want to plot the responses, as in for a paper
    fig, ax = plt.subplots(1, 2, sharex=true)
    plotcurve(V0, ϕ, q, fig, gridalpha=gridalpha)
    fig
end   

function plotcurve(V0, ϕ, q, fig; iaxis=1, gridalpha=0.5, modelalpha=0.5, 
    lcolor="nocolor", writelabel=true)
    # plotting responses into an existing figure, useful for plotting multiple responses
    ax = fig.axes
    if lcolor == "nocolor" # as in use default PyPlot colors
        ax[iaxis].plot(q, V0)
        ax[iaxis+1].plot(q, ϕ)
    else # plot with specified color
        ax[iaxis].plot(q, V0, color=lcolor, alpha=modelalpha)
        ax[iaxis+1].plot(q, rem2pi.(ϕ, RoundNearest), color=lcolor, alpha=modelalpha)
    end
    if writelabel # useful on the last iteration of model plotting
        writelabels(ax, iaxis, gridalpha)
        ax[iaxis+1].set_ylim(-pi, pi)
    end    
    fig.tight_layout()
    nothing
end  

function writelabels(ax, iaxis, gridalpha)
    ax[iaxis].set_xlabel("Pulse moment A-s")
    ax[iaxis].set_ylabel("Amplitude")
    ax[iaxis].set_xscale("log")
    ax[iaxis].grid(b=true, which="both", alpha=gridalpha)
    ax[iaxis+1].set_xlabel("Pulse moment A-s")
    ax[iaxis+1].set_ylabel("phase")
    ax[iaxis+1].grid(b=true, which="both", alpha=gridalpha)
end  

function plot_model_field(S::SMRSounding, opt::Options, optn::OptionsNuisance; 
        gridalpha=0.5, figsize=(10,4), lcolor="nocolor", modelalpha=0.5, burninfrac=0.5, 
        decfactor=10)
    fig = figure(figsize=(figsize))
    s1 = subplot(131)
    s2 = subplot(132)
    s3 = subplot(133, sharex=s2)
    M = assembleTat1(opt, :fstar, temperaturenum=1, burninfrac=burninfrac)[1:decfactor:end]
    Mnu = assemblenuisancesatT(optn, burninfrac = burninfrac, temperaturenum = 1)[1:decfactor:end]
    T = S.linearsat ? x->x : x->10^x
    for (imodel, (m, nu)) in enumerate(zip(M, Mnu))
        w = T.(m)[:]
        response = SNMRForward.forward(S.Fm, w)
        Vresp = abs.(response)
        ϕresp = angle.(response) .+ nu[1]
        writelabel = imodel == length(M) ? true : false
        plotmodelcurve(w, S.Fm.zgrid, Vresp, ϕresp, S.Fm.qgrid, fig, writelabel=writelabel,
            gridalpha=gridalpha, lcolor=lcolor, modelalpha=modelalpha)
    end
    plotdata(S, fig, iaxis=2, gridalpha=gridalpha)    
end  

end