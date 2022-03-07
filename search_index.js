var documenterSearchIndex = {"docs":
[{"location":"forward_modelling/#Forward-Modelling","page":"Forward Modelling","title":"Forward Modelling","text":"","category":"section"},{"location":"forward_modelling/","page":"Forward Modelling","title":"Forward Modelling","text":"SNMRForward defines functions and objects to be used for forward modelling of surface NMR free-induction decay (FID) experiments.","category":"page"},{"location":"forward_modelling/","page":"Forward Modelling","title":"Forward Modelling","text":"MRSForward\n\nMRSForward_square\n\nforward","category":"page"},{"location":"forward_modelling/#SNMRForward.MRSForward","page":"Forward Modelling","title":"SNMRForward.MRSForward","text":"MRSForward(R::Real, zgrid::AbstractVector{<:Real},\nqgrid::AbstractVector{<:Real}, ϕ::Real, Be::Real,\ncondLEM::ConductivityModel; nrvals = 200, temp=300.0,\nqwe=true)\n\nReturns an MRSForward forward modelling object, with the kernel computed for a circular loop.\n\nParameters:\n\nR: loop radius in metres\nzgrid: depth cell boundaries, in metres below surface.\nqgrid: pulse moments used for the experiment, in Ampere-seconds.\nϕ: magnetic field inclination in radians.\nBe: Earth field magnitude in Tesla.\ncondLEM: layered-earth conductivity model, as a ConductivityModel object.\n\nOptional parameters:\n\nnrvals: number of radial points to use to evaluate the horizontal integration of the kernel.\ntemp: temperature in Kelvin, affects the magnitude of equilibrium magnetisation.\nqwe: whether to use quadrature with extrapolation (QWE) to evaluate Hankel transforms for the kernel. This is slower, but more accurate at shallower depths, than the alternative 801-point digital filter.\n\n\n\n\n\n","category":"type"},{"location":"forward_modelling/#SNMRForward.MRSForward_square","page":"Forward Modelling","title":"SNMRForward.MRSForward_square","text":"MRSForward_square(L::Real, zgrid::AbstractVector{<:Real},\nqgrid::AbstractVector{<:Real}, ϕ::Real, θ::Real, Be::Real,\ncondLEM::ConductivityModel; nxpoints = 128, temp=300.0)\n\nReturns an MRSForward forward modelling object, with the kernel computed for a square loop.\n\nParameters:\n\nL: loop radius in metres\nzgrid: depth cell boundaries, in metres below surface.\nqgrid: pulse moments used for the experiment, in Ampere-seconds.\nϕ: magnetic field inclination in radians.\nθ: horizontal angle between loop orientation and magnetic north, in radians.\nBe: Earth field magnitude in Tesla.\ncondLEM: layered-earth conductivity model, as a ConductivityModel object.\n\nOptional parameters:\n\nnxpoints: number of points to use in each dimension for the horizontal integration of the kernel.\ntemp: temperature in Kelvin, affects the magnitude of equilibrium magnetisation.\n\n\n\n\n\n","category":"function"},{"location":"forward_modelling/#SNMRForward.forward","page":"Forward Modelling","title":"SNMRForward.forward","text":"forward(F::MRSForward, m::Vector{<:Real})\n\nComputes the complex response amplitude of a water content model for a given FID experiment setup, described by an MRSForward object. Returns response as a complex vector, in units of volts, one  element per pulse moment.\n\nParameters:\n\nF: an MRSForward object containing the kernel for modelling\nm: a vector containing water saturation at each depth in F.zgrid.\n\n\n\n\n\n","category":"function"},{"location":"#SNMRForward.jl-Documentation","page":"SNMRForward.jl Documentation","title":"SNMRForward.jl Documentation","text":"","category":"section"},{"location":"","page":"SNMRForward.jl Documentation","title":"SNMRForward.jl Documentation","text":"","category":"page"}]
}
