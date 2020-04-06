
@testset "symmetry equivalent wavenumbers" begin

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# ms = MultipleScattering # just in case Circle and Sphere conflict with definitions from other packages.

s1 = Specie(
    Acoustic(spatial_dim; ρ=0.2, c=0.1), Sphere(0.4);
    volume_fraction=0.2
)
species = [s1]

ω = 0.9

AP_det = dispersion_equation(ω, medium, species, PlanarAzimuthalSymmetry())
P_det = dispersion_equation(ω, medium, species, PlanarSymmetry())
AR_det = dispersion_equation(ω, medium, species, AzimuthalSymmetry())
R_det = dispersion_equation(ω, medium, species, WithoutSymmetry())

AP_kps = wavenumbers(ω, medium, species;
    num_wavenumbers = 6, mesh_points = 4,
    symmetry = PlanarAzimuthalSymmetry())

P_kps = wavenumbers(ω, medium, species;
    num_wavenumbers = 6, mesh_points = 4, k_effs = AP_kps,
    symmetry = PlanarSymmetry())

# As plane waves with azimuthal symmetry is a sub-case of plane-waves, and all materials allow for the effective wavenumbers of plane waves, all the below determinant equations should be satisfied
@test maximum(AP_det.(AP_kps)) < 1e-4
@test maximum(P_det.(AP_kps)) < 1e-4
@test maximum(AR_det.(AP_kps)) < 1e-10
@test maximum(R_det.(AP_kps)) < 1e-20

# However, there do exist effective wavenumbers for plane-waves which have eigen-vectors that do not satisfy azimuthal symmetry. This is why maximum(AP_det.(P_kps)) != 0.0
@test maximum(P_det.(P_kps)) < 1e-4
@test maximum(AR_det.(P_kps)) < 1e-10
@test maximum(R_det.(P_kps)) < 1e-20


## Test low volume fraction

volfrac = 0.02
s1 = Specie(
    Acoustic(spatial_dim; ρ=2.2, c=3.1), Sphere(0.4);
    volume_fraction = volfrac / 2.0
)

s2 = Specie(
    Acoustic(spatial_dim; ρ=12.2, c=23.1), Sphere(0.2);
    volume_fraction = volfrac / 2.0
)
species = [s1,s2]

ω = 1.9

kp_lowvol = wavenumber_low_volumefraction(ω, medium, species; basis_order=2)

AP_kps = wavenumbers(ω, medium, species; basis_order=2,
    num_wavenumbers = 2, mesh_size = 2.0, mesh_points = 4,
    symmetry = PlanarAzimuthalSymmetry())

@test abs(AP_kps[1] - kp_lowvol) / abs(kp_lowvol) < 10 * volfrac^3


# Choose specific materials to calculate eignvectors

# θin = 0.0
# basis_order = 1
#
# normal = [-1.0,0.0] # an outward normal to the surface
# material = Material(ms.Sphere(4.0),species)
# source = PlaneSource(medium, [cos(θin),0.0,sin(θin)])
#
# MM_regular = eigensystem(ω, source, material;
#         basis_order = basis_order,
#         basis_order_field = 2*basis_order,
# )


end
