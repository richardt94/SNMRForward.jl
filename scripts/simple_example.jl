using SpecialFunctions, Hankel, Elliptic, PyPlot

#constants in SI units
mu_0 = 4*pi*10^-7

#frequency-domain field in horizontal wavenumber space

I = 1; #current (amps)
R = 10; #radius (metres)

Hz_wavenumber(κ,z) = I * R * besselj1(κ*R) * exp(- κ * abs(z))/2;

##
max_radius = 3*R;
hank_0 = QDHT(0,max_radius,201);

##
Hz_radius = hank_0 \ Hz_wavenumber.(hank_0.k, R)

##
zgrid = 0:0.01*R:10*R
Hzk_mat = [Hz_wavenumber(κ,z) for κ in hank_0.k, z in zgrid]
Hzr_mat = hank_0 \ Hzk_mat

plt.figure();
plt.pcolor(hank_0.r, zgrid, Hzr_mat', vmin = -0.06, vmax = 0.06);
gca().invert_yaxis();
xlabel("r (m)")
ylabel("z (m)")
display(gcf());
plt.close("all");

##

#analytic method for spatial fields
#(without numerical Hankel transform)
function Hz_analytic(r,z)
    αsq = R^2 + r^2 + z^2 - 2*r*R;
    βsq = R^2 + r^2 + z^2 + 2*r*R;
    β = sqrt(βsq);
    ksq = 1 - αsq/βsq;
    C = 1/pi;
    
    C/(2*αsq*β)*((R^2 - r^2 - z^2)*E(ksq) + αsq*K(ksq))
end

Hz_analytic_mat = [Hz_analytic(r,z) for r in hank_0.r, z in zgrid]

##
plt.figure();
plt.pcolor(hank_0.r, zgrid, Hz_analytic_mat', vmin = -0.06, vmax = 0.06);
gca().invert_yaxis();
xlabel("r (m)")
ylabel("z (m)")
display(gcf());
plt.close("all");

#evaluate "B-response" for each layer using the conductivities
#and at spatial wavenumbers κ
##

RealVector = Vector{T} where T <: Real
function responses(σ::RealVector, d::RealVector, κ::RealVector, ωl::Real)
    B = zeros(ComplexF64, length(κ), length(σ));
    α = [sqrt(k^2 - im * ωl * mu_0 * σm)  for k in κ, σm in σ]
    B[:,end] = -α[:,end]
    m = size(B,2) - 1
    print(m)
    while m > 0
        thamdm = tanh.(α[:,m]*d[m])
        B[:,m] = α[:,m] .* (B[:,m+1] .- α[:,m].*thamdm)./(α[:,m] - B[:,m+1].*thamdm)
        m -= 1
    end
    return B, α
end


##
dipole_approximation = [25/((r^2+z^2)^(3/2)) * (3*z^2/(z^2 + r^2) - 1) for r in hank_0.r, z in zgrid]

##
plt.figure();
plt.pcolor(hank_0.r, zgrid, dipole_approximation', vmin = -0.06, vmax = 0.06);
gca().invert_yaxis();
xlabel("r (m)")
ylabel("z (m)")
display(gcf());
plt.close("all");

##
#far-field comparison
figure()
plot(zgrid[end-200:end], dipole_approximation[1,end-200:end])
plot(zgrid[end-200:end], Hzr_mat[1,end-200:end])
plot(zgrid[end-200:end], Hz_analytic_mat[1,end-200:end])
legend(["dipole", "numerical hankel transform", "analytic"])
xlabel("z (m)")
ylabel("Hz field (1/m)")
display(gcf())
close("all")

##
ωl = 2.0e3 #Hz, typical for Earth's field strength
B_halfspace, alpha_halfspace = responses([0,0.02], [Inf,Inf], hank_0.k, ωl)

#Schelkunoff potential at each z
function phiz(phi0, Bresponse, α, d, zgrid)
    #schelkunoff potential at the top of each layer
    phi_tops = zeros(ComplexF64, size(Bresponse)...)
    phi_tops[:,1] = phi0
    for m=1:size(phi_tops,2)-1
        phi_tops[:,m+1] = phi_tops[:,m] .* (α[:,m] + Bresponse[:,m])./(α[:,m] + Bresponse[:,m+1]) .* exp.(α[:,m]*d[m])
    end
    #convert depth to z-values
    z_interface = cumsum(d)

    phiz = zeros(ComplexF64, length(phi0), length(zgrid))
    for z in zgrid
        
    end
    

end