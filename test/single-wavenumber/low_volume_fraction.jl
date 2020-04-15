
@testset "low volume fraction 3D acoustics" begin

## Test low volume fraction

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

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

end
