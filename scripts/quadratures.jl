using FastGaussQuadrature, SpecialFunctions, PyPlot, Elliptic
##
x, w = gausslegendre(8)
f(x) = x^4
transpose(w)*f.(x) ≈ 2/5
##
function quad_int(f, a, b)
    u = (b - a)/2 * (x .+ 1) .+ a
    (b-a)/2 * transpose(w) * f.(u)
end
##
quad_int(x->x^2, 2, 3)
##
Hz_wavenumber(k) = besselj1(k)/2
##
b_roots = approx_besselroots(1,300)
##
function psum_kern(fk, kern, r, R)
    order = length(b_roots)
    psums = zeros(order)

    integrand(k) = fk(k) * k * kern(k*r)
    
    a = 0
    b = b_roots[1]/R
    psums = NaN * ones(order)
    psums[1] = quad_int(integrand, a, b)
    epsilon = NaN * ones(order, order)
    function calc_eps(i,j)
        ~isnan(epsilon[i+1,j+1]) && return epsilon[i+1,j+1]
        if j == 0
            if isnan(psums[i+1])
                a = b
                b = b_roots[i+1]/R
                psums[i+1] = psums[i] + quad_int(integrand, a, b)
            end
            epsilon[i+1,1] = psums[i+1]
        elseif j == 1
            epsilon[i+1,2] = 1/(calc_eps(i+1,0) - calc_eps(i,0))
        else
            epsilon[i+1,j+1] = calc_eps(i+1, j-2) + 1/(calc_eps(i+1,j-1) - calc_eps(i,j-1))
        end
        return epsilon[i+1,j+1]
    end

    os = 0
    s = 0
    tol = 1e-6
    for i=0:(order-1)
        if i%2 == 0
            s = calc_eps(0,i)
        else
            s = calc_eps(1,i-1)
        end

        abs((s - os)/s) < tol && break
        
        os = s
    end
    s
end
##
rgrid = 0.01:0.01:2
# figure()
# plot(rgrid, psum_j0.(Hz_wavenumber, rgrid))
##
# gcf()
##
# close("all")
##
result = psum_kern.(Hz_wavenumber, besselj0, rgrid, 1)
##
function Hz_analytic(r)
    ksq = 1 - ((r-1)/(r+1))^2
    ((r-1)*K(ksq) - (r+1)*E(ksq))/(2*pi*(r+1)*(r-1))
end
##
res_true = Hz_analytic.(rgrid)
##
figure()
plot(rgrid, result, rgrid, res_true)
ylim([-20,20])
##
gcf()
##
function Hz_analytic(r,z)
    αsq = 1 + r^2 + z^2 - 2*r;
    βsq = 1 + r^2 + z^2 + 2*r;
    β = sqrt(βsq);
    ksq = 1 - αsq/βsq;
    C = 1/pi;
    
    C/(2*αsq*β)*((1 - r^2 - z^2)*E(ksq) + αsq*K(ksq))
end

function Hr_analytic(r,z)
    αsq = 1 + r^2 + z^2 - 2*r;
    βsq = 1 + r^2 + z^2 + 2*r;
    β = sqrt(βsq);
    ksq = 1 - αsq/βsq;
    C = 1/pi;
    
    C * z / (2*αsq*β*r) * ((1 + r^2 + z^2)*E(ksq) - αsq*K(ksq))
end
##

Hz_wavenumber(k,z) = exp(-k*z) * Hz_wavenumber(k)
Hr_wavenumber(k) = besselj1(k)/2
Hr_wavenumber(k,z) = exp(-k*z) * Hz_wavenumber(k)

##
result_Hr = psum_kern.(Hr_wavenumber, besselj1, rgrid,1)
result_analytic_Hr = Hr_analytic.(rgrid,0)
##
figure()
plot(rgrid, result_Hr)
plot(rgrid, result_analytic_Hr)
gcf()
## 2D results
zgrid = 10 .^(-2:0.05:0.5)
Hz2d = reduce(hcat, collect(psum_kern.(k -> Hz_wavenumber(k,z), besselj0, rgrid, 1) for z in zgrid))
Hr2d = reduce(hcat, collect(psum_kern.(k -> Hr_wavenumber(k,z), besselj1, rgrid, 1) for z in zgrid))
Hza2d = [Hz_analytic(r,z) for r in rgrid, z in zgrid]
Hra2d = [Hr_analytic(r,z) for r in rgrid, z in zgrid]
##
clip = 0.4
fig,ax = subplots(2,2, figsize=(10,10))
sca(ax[1,1])
pcolor(rgrid, zgrid, transpose(Hz2d), vmin=-clip, vmax=clip)
gca().invert_yaxis()
title("QWE Hz")
sca(ax[1,2])
pcolor(rgrid, zgrid, transpose(Hza2d), vmin=-clip, vmax=clip)
gca().invert_yaxis()
title("analytic Hz")
sca(ax[2,1])
pcolor(rgrid, zgrid, transpose(Hr2d), vmin=-clip, vmax=clip)
gca().invert_yaxis()
title("QWE Hr")
sca(ax[2,2])
pcolor(rgrid, zgrid, transpose(Hra2d), vmin=-clip, vmax=clip)
gca().invert_yaxis()
title("analytic Hr")
for a in ax[:]
    sca(a)
    xlabel("r/R")
    ylabel("z/R")
end
gcf()
##
# design a way of actually evaluating a messier forward with QWE
# starting with params needed for LE formulation
using SNMRForward
d = Vector{Float64}()
σ = [0.001]
ωl = 2*π*2.5e3
R = 50
rgrid = 0.01*R:0.02*R:4*R
zgrid = R * 10 .^(-2:0.05:0.4)
##
function calc_H(rs, zs, R, σ, d, ωl; quad_order = 8, max_zero = 50, tol=1e-7)
    x, w = gausslegendre(quad_order)
    b_roots = [0; approx_besselroots(1, 150)/R]
    scaler = (b_roots[2:end] .- b_roots[1:end-1])/2
    nr = length(rs)
    nz = length(zs)
    Hz = zeros(Complex, nr,nz)
    Hr = zeros(Complex, nr,nz)
    ks = reduce(vcat, [(x .+ 1)/2 * (b-a) .+ a for (a,b) in zip(b_roots[1:end-1], b_roots[2:end])])
    
    #precompute layer responses for each k
    B, α = SNMRForward.responses(σ, d, ks, ωl)
    phi, phip = SNMRForward.phi_coeffs(B, α, d, zs)
    #reflection coefficient at first interface
    rte = 2 * ks ./ (ks .+ B[:,1])
    ##
    for (ir, r) = enumerate(rs)
        #integrand to compute magnetic field in free space at z = 0
        free_integrand_z(k) = k * SNMRForward.H_free(k, R) * besselj0(k*r)
        free_integrand_r(k) = - SNMRForward.H_free(k, R) * besselj1(k*r)
        for iz = 1:nz
            psums = NaN * ones(Complex, max_zero, 2)
            epsilon = NaN * ones(Complex, max_zero, max_zero, 2)
            function integrate(iq, isz)
                slice = (iq*quad_order+1):(iq+1)*quad_order
                scaler[iq+1] * transpose(w .* rte[slice]) * (isz ?
                    (free_integrand_z.(ks[slice]) .* phi[slice,iz])
                    :
                    (free_integrand_r.(ks[slice]) .* phip[slice,iz]))
            end
            psums[1,1] = integrate(0, false)
            psums[1,2] = integrate(0, true)

            function calc_eps(i,j,isz)
                ~isnan(epsilon[i+1,j+1,isz+1]) && return epsilon[i+1,j+1,isz+1]
                if j == 0
                    if isnan(psums[i+1,isz+1])
                        psums[i+1,isz+1] = psums[i,isz+1] + integrate(i,isz)
                    end
                    epsilon[i+1,1,isz+1] = psums[i+1,isz+1]
                elseif j == 1
                    epsilon[i+1,2,isz+1] = 1/(calc_eps(i+1,0,isz) - calc_eps(i,0,isz))
                else
                    epsilon[i+1,j+1,isz+1] = 
                        calc_eps(i+1,j-2,isz) + 1/(calc_eps(i+1,j-1,isz) - calc_eps(i,j-1,isz))
                end
                return epsilon[i+1,j+1,isz+1]
            end

            function adaptive_shanks(isz)
                os = 0
                s = 0
                for i=0:(max_zero-1)
                    if i%2 == 0
                        s = calc_eps(0,i,isz)
                    else
                        s = calc_eps(1,i-1,isz)
                    end
            
                    abs((s - os)/s) < tol && break
                    
                    os = s
                end
                s
            end

            Hz[ir,iz] = adaptive_shanks(true)
            Hr[ir,iz] = adaptive_shanks(false)
        end
    end

    Hz, Hr
end

##
@time Hz, Hr = calc_H(rgrid, zgrid, R, σ, d, ωl; quad_order=16)
##
@time (Hz_df, Hr_df) = SNMRForward.magfields(R, ωl, σ, d, rgrid, zgrid)
##
Hz_rd = abs.((Hz .- Hz_df)./Hz)
Hr_rd = abs.((Hr .- Hr_df)./Hr)
##
fig, ax = subplots(3,2, figsize=(10,15))
sca(ax[1,1])
pcolor(rgrid, zgrid, real.(transpose(Hz)), vmin=-0.02, vmax=0.02)
title("Hz QWE")
sca(ax[1,2])
pcolor(rgrid, zgrid, real.(transpose(Hz_df)), vmin=-0.02, vmax=0.02)
title("Hz dig filter")
sca(ax[2,1])
pcolor(rgrid, zgrid, real.(transpose(Hr)), vmin=-0.02, vmax=0.02)
title("Hr QWE")
sca(ax[2,2])
pcolor(rgrid, zgrid, real.(transpose(Hr_df)), vmin=-0.02, vmax=0.02)
title("Hr dig filter")
sca(ax[3,1])
pcolor(rgrid, zgrid, real.(transpose(Hz_rd)), vmin=0, vmax=0.05)
title("Hz difference")
sca(ax[3,2])
pcolor(rgrid, zgrid, real.(transpose(Hr_rd)), vmin=0, vmax=0.05)
title("Hr difference")


for a in ax[:]
    sca(a)
    a.invert_yaxis()
    xlabel("r (m)")
    ylabel("z (m)")
end
gcf()
##
# 1d kernels
qs = [0.1,0.25,0.5,0.75,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
k1d = reduce(hcat, [SNMRForward.kernel_1d(q, 65*π/180, ωl, Hz, Hr, rgrid; n_theta_points=250) for q in qs])
k1d_df = reduce(hcat, [SNMRForward.kernel_1d(q, 65*π/180, ωl, Hz_df, Hr_df, rgrid; n_theta_points=250) for q in qs])
##
figure()
plot(k1d[:,5],zgrid)
plot(k1d_df[:,5],zgrid)
ylim([1,maximum(zgrid)])
gca().invert_yaxis()
gca().set_yscale("log")
gcf()
gcf()
##
close("all")
##
# forward modelling with qwe and digital filter kernels
m0 = SNMRForward.mag_factor(300) * ωl / SNMRForward.γh
F = SNMRForward.MRSForward(m0 * k1d, qs, zgrid, zgrid[2:end] .- zgrid[1:end-1])
F_df = SNMRForward.MRSForward(m0 * k1d_df, qs, zgrid, zgrid[2:end] .- zgrid[1:end-1])
##
w = zeros(length(zgrid))
w[(zgrid .>= 10) .& (zgrid .<= 20)] .= 1
w
V = SNMRForward.forward(F,w)
V_df = SNMRForward.forward(F_df, w)
##
figure()
plot(qs, real.(V))
plot(qs, real.(V_df))
gcf()
##