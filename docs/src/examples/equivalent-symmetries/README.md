# Equivalent symmetries

Here we show how to formulate dispersion equations for wave mode which are not just plane waves. See [Background](@ref) for an overview.

The dispersion equations become simpler the more symmetries that are shared between the incident wave and the material geometry. The simplest case is [`PlanarAzimuthalSymmetry`](@ref), and includes the case where a plane wave is directly incident upon a flat surface. If no symmetries are present, then the type [`WithoutSymmetry`](@ref) is used and leads to a general dispersion equation for materials which occupy a [simple connected domain](https://en.wikipedia.org/wiki/Simply_connected_space). This general dispersion equation has spurious roots and is computationally heavier to solve, see [Gower & Kristensson 2020](https://arxiv.org/pdf/2010.00934.pdf).  

## Choose the microstructure
```julia
spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# Below we explicitly define the shape of the particles as being spheres
s1 = Specie(
    Acoustic(spatial_dim; ρ=1.0, c=0.5), Sphere(spatial_dim, 0.4);
    volume_fraction=0.2
);

# We use just one specie to speed up the calculations.
species = [s1]
```

## Calculate the wavenumbers
```julia
ω = 0.9
tol = 1e-7

# Calculate the wavenumbers for PlanarAzimuthalSymmetry()
AP_kps = wavenumbers(ω, medium, species;
    num_wavenumbers = 4, tol = tol,
    symmetry = PlanarAzimuthalSymmetry())

# Calculate the wavenumbers for just PlanarSymmetry()
P_kps = wavenumbers(ω, medium, species;
    num_wavenumbers = 4, tol = tol,
    symmetry = PlanarSymmetry{3}())

# Select a subset to test. The dispersion equation WithoutSymmetry can be unstable for higher order effective wavenumbers.
AP_kps = AP_kps[1:min(length(AP_kps),14)]
P_kps = P_kps[1:min(length(P_kps),14)]
```

## Test different dispersion equations
```julia
AP_det = dispersion_equation(ω, medium, species, PlanarAzimuthalSymmetry{spatial_dim}())
P_det = dispersion_equation(ω, medium, species, PlanarSymmetry{spatial_dim}())
AR_det = dispersion_equation(ω, medium, species, AzimuthalSymmetry{spatial_dim}())
R_det = dispersion_equation(ω, medium, species, WithoutSymmetry{spatial_dim}())

# As plane waves with azimuthal symmetry is a sub-case of plane-waves, and all materials allow for the effective wavenumbers of plane waves, all the below determinant equations should be satisfied
maximum(AP_det.(AP_kps)) < tol
maximum(P_det.(AP_kps)) < tol
maximum(AR_det.(AP_kps)) < tol^2
maximum(R_det.(AP_kps)) < tol^3

# However, there do exist effective wavenumbers for plane-waves which have eigen-vectors that do not satisfy azimuthal symmetry. This is why maximum(AP_det.(P_kps)) != 0.0

maximum(AP_det.(P_kps)) > tol
maximum(P_det.(P_kps)) < tol
maximum(AR_det.(P_kps)) < tol^2
maximum(R_det.(P_kps)) < tol^3
```
