using EffectiveWaves, Test, LinearAlgebra

# 2.2159918607732766e-10
#  3.6428008939057475e-8
#  2.1279008239602572e-7


@testset "low frequency 3D acoustics" begin

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# ms = MultipleScattering # just in case Circle and Sphere conflict with definitions from other packages.

s1 = Specie(
    Acoustic(spatial_dim; ρ=10.2, c=10.1), Sphere(0.001);
    volume_fraction=0.2
);
s2 = Specie(
    Acoustic(spatial_dim; ρ=0.2, c=4.1), Sphere(0.0002);
    volume_fraction=0.15
);
species = [s1,s2]

ωs = LinRange(1e-2,1e-1,3)
tol = 1e-7
basis_order = 1

θ = 0.0
psource = PlaneSource(medium, [sin(θ),0.0,cos(θ)]);
source = plane_source(medium; direction = [sin(θ),0.0,cos(θ)])
material = Material(Sphere(4.0),species);

basis_field_order = 3

ks = ωs ./ medium.c

opts = Dict(
    :tol => tol, :num_wavenumbers => 2,
    :mesh_size => 2.0, :mesh_points => 10,
    :basis_order => basis_order, :basis_field_order => basis_field_order
);
AP_kps = [
    wavenumbers(ω, medium, species; symmetry = PlanarAzimuthalSymmetry(), opts...)
for ω in ωs]

k_effs = [kps[1] for kps in AP_kps]

eff_medium = effective_medium(medium, species)

k_lows = ωs ./ eff_medium.c

errs = abs.(k_effs - k_lows)
@test errs[1] < errs[2] < errs[3]
@test maximum(errs) < tol
@test errs[1] < 10.0 * tol^2

A_waves = [
    WaveMode(ωs[i], k_effs[i], psource, material;
        basis_order = basis_order, basis_field_order = basis_field_order)
for i in eachindex(k_effs)]

scat_azis = material_scattering_coefficients.(A_waves);

r = maximum(outer_radius.(species))
material_low = Material(Sphere(4.0 - r),species);
effective_sphere = Particle(eff_medium, material_low.shape)

Linc = basis_field_order + basis_order;
n_to_l = [l for l = 0:Linc for m = -l:l];

errs =  map(eachindex(ωs)) do i
    source_coefficients =  regular_spherical_coefficients(source)(Linc,zeros(3),ωs[i])
    Tmat = MultipleScattering.t_matrix(effective_sphere, medium, ωs[i], Linc)
    norm(scat_azis[i] - diag(Tmat)[n_to_l .+ 1] .* source_coefficients) / norm(source_coefficients)
end

@test sum(errs .< maximum(outer_radius.(species)) .* abs.(ks)) == length(ks)
@test errs[1] < tol
@test maximum(errs) < 5.0 * tol

end
