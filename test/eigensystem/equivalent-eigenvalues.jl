using EffectiveWaves, Test, LinearAlgebra

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

ms = MultipleScattering # just in case Circle and Sphere conflict with definitions from other packages.

s1 = Specie(
    Acoustic(spatial_dim; ρ=0.2, c=0.1),ms.Sphere(0.4);
    volume_fraction=0.2
)
species = [s1]

ω = 0.8
θin = 0.0
basis_order = 1

normal = [-1.0,0.0] # an outward normal to the surface
material = Material(ms.Sphere(4.0),species)
source = PlaneSource(medium, [cos(θin),0.0,sin(θin)])

MM_regular = eigensystem(ω, source, material;
        basis_order = basis_order,
        basis_order_field = 2*basis_order,
)
