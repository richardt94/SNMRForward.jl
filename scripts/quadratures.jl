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
close("all")
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
zgrid = 10 .^(-2:0.05:0.5)
##
function calc_H(rs, zs, R, σ, d, ωl; quad_order = 8, max_zero = 150)
    x, w = gausslegendre(quad_order)
    b_roots = [0; approx_besselroots(1, 150)/R]
    scaler = (b_roots[2:end] .- b_roots[1:end-1])/2
    nr = length(rs)
    nz = length(zs)
    Hz = zeros(Complex, nr,nz)
    ks = reduce(vcat, [(x .+ 1)/2 * (b-a) .+ b for (a,b) in zip(b_roots[1:end-1], b_roots[2:end])])
    #precompute layer responses for each k
    B, α = SNMRForward.responses(σ, d, ks, ωl)
    phi, phip = SNMRForward.phi_coeffs(B, α, d, zs)
    ##
    for (ir, r) = enumerate(rs)
        #integrand to compute magnetic field in free space at z = 0
        free_integrand(k) = SNMRForward.H_free(k, R) * besselj0(k*r)
        for (iz, z) = enumerate(zs)
            psums = NaN * ones(Complex, max_zero)
            epsilon = NaN * ones(Complex, max_zero, max_zero)
            function integrate(iq)
                slice = (iq*quad_order+1):(iq+1)*quad_order
                scaler[iq+1] * transpose(w) * 
                    (free_integrand.(ks[slice]) .* phi[slice,iz])
            end
            psums[1] = integrate(0)

            function calc_eps(i,j)
                ~isnan(epsilon[i+1,j+1]) && return epsilon[i+1,j+1]
                if j == 0
                    if isnan(psums[i+1])
                        psums[i+1] = psums[i] + integrate(i)
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
            for i=0:(max_zero-1)
                if i%2 == 0
                    s = calc_eps(0,i)
                else
                    s = calc_eps(1,i-1)
                end
        
                abs((s - os)/s) < tol && break
                
                os = s
            end
            Hz[ir,iz] = s

        end
    end

    Hz    
end

##
Hz = calc_H(0.01:0.02:2, zgrid, R, σ, d, ωl)
##
