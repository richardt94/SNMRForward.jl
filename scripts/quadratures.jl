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
function psum_j0(fk, r)
    order = length(b_roots)
    psums = zeros(order)

    integrand(k) = fk(k) * k * besselj0(k*r)
    a = 0
    b = b_roots[1]
    psums[1] = quad_int(integrand, a, b)
    for i in 2:order
        a = b
        b = b_roots[i]
        psums[i] = psums[i-1] + quad_int(integrand, a, b)
    end

    epsilon = NaN * ones(order, order)
    function calc_eps(i,j)
        ~isnan(epsilon[i+1,j+1]) && return epsilon[i+1,j+1]
        if j == 0
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
rgrid = 0.01:0.01:4
figure()
plot(rgrid, psum_j0.(Hz_wavenumber, rgrid))
##
gcf()
##
close("all")
##
result = psum_j0.(Hz_wavenumber, rgrid)
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