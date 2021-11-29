#kernel function evaluated at a point, given co- and counter-rotating fields
#and phase lag
#must be scaled by equilibrium magnetisation (depends on magnitude of earth field)
point_kernel(q, Bplus, Bminus, ζ, ωl) = 
    2 * ωl * sin(γh * q * Bplus) * Bminus * exp(-2*im*ζ)

function kernel_1d(q, ϕ, ωl, Hz, Hr, rgrid; n_theta_points = 100)
    #this uses trapezoidal integration in theta and r
    thetagrid = range(0, 2*pi, length=n_theta_points)
    #radial integral scale
    dr = rgrid[2:end] .- rgrid[1:end-1]
    #azimuthal integral scale
    dtheta = thetagrid[2] - thetagrid[1]

    k1d = zeros(ComplexF64, size(Hz,2))

    for (i_th, θ) = enumerate(thetagrid)
        Hparams = SNMRForward.co_counter_field.(Hz, Hr, ϕ, θ)

        Hco = first.(Hparams)
        ζ = last.(Hparams)
        Hctr = reshape([a[2] for a in Hparams[:]], size(Hco)...)

        kernel = SNMRForward.point_kernel.(q, μ0 * Hco, μ0 * Hctr, ζ, ωl)
        rkr = rgrid .* kernel
        if i_th == 1 || i_th == length(thetagrid)
            k1d += 0.25 * dtheta*transpose(rkr[1:end-1,:] .+ rkr[2:end,:])*dr
        else
            k1d += 0.5 * dtheta*transpose(rkr[1:end-1,:] .+ rkr[2:end,:])*dr
        end
    end
    k1d
end