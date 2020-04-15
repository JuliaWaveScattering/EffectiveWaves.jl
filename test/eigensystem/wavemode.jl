using EffectiveWaves, Test, LinearAlgebra

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# ms = MultipleScattering # just in case Circle and Sphere conflict with definitions from other packages.

s1 = Specie(
    Acoustic(spatial_dim; ρ=0.2, c=0.1), Sphere(0.4);
    volume_fraction=0.2
);
species = [s1]

ω = 0.9
tol = 1e-6

AP_det = dispersion_equation(ω, medium, species, PlanarAzimuthalSymmetry())
P_det = dispersion_equation(ω, medium, species, PlanarSymmetry())
AR_det = dispersion_equation(ω, medium, species, AzimuthalSymmetry())
R_det = dispersion_equation(ω, medium, species, WithoutSymmetry())

AP_kps = wavenumbers(ω, medium, species;
    num_wavenumbers = 6, tol = tol,
    symmetry = PlanarAzimuthalSymmetry())

θ = 0.0
normal = [0.0,0.0,-1.0]; # an outward normal to the surface
material = Material(Sphere(4.0),species);
source = PlaneSource(medium, [sin(θ),0.0,cos(θ)]);

k_eff = AP_kps[1]

MM = eigensystem(ω, source, material; basis_order = 2)

MM_svd = svd(MM(k_eff))

inds = findall(MM_svd.S .< tol)

inds = findall(MM_svd.S .< -1)
