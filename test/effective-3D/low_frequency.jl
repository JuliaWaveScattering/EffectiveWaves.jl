using EffectiveWaves, Test, LinearAlgebra

@testset "low frequency 3D acoustics" begin

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# NOTE: we so far only understand the region of the low frequency effective medium when all the particles have the same radius. That is, we get exact matches below (Monte Carlo also confirms this). However, when using species with different radius, it is not clear what effective bounding shape to use.

s1 = Specie(
    Acoustic(spatial_dim; ρ=10.2, c=10.1), Sphere(0.001);
    volume_fraction=0.2
);
s2 = Specie(
    Acoustic(spatial_dim; ρ=0.2, c=4.1), Sphere(0.001);
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

basis_field_order = 6

ks = ωs ./ medium.c

opts = Dict(
    :tol => tol, :num_wavenumbers => 2,
    :mesh_size => 2.0, :mesh_points => 10,
    :basis_order => basis_order, :basis_field_order => basis_field_order
);
AP_kps = [
    wavenumbers(ω, medium, species;
        symmetry = PlanarAzimuthalSymmetry(),
        numberofparticles = material.numberofparticles,
        opts...
    )
for ω in ωs]

k_effs = [kps[1] for kps in AP_kps]

eff_medium = effective_medium(medium, species; numberofparticles = material.numberofparticles)

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
material_low = Material(Sphere(outer_radius(material.shape) - r),species);
effective_sphere = Particle(eff_medium, material_low.shape)

Linc = basis_field_order + basis_order;

# I can not work out theoretically why this scale_number_density should be needed below.
scale_number_density = 1.0 - 1.0 / material.numberofparticles

errs =  map(eachindex(ωs)) do i
    source_coefficients =  regular_spherical_coefficients(source)(Linc,zeros(3),ωs[i])
    Tmat = MultipleScattering.t_matrix(effective_sphere, medium, ωs[i], Linc)
    Tmat = Tmat ./ scale_number_density
    norm(scat_azis[i] - Tmat * source_coefficients) / norm(scat_azis[i])
end

@test sum(errs .< maximum(outer_radius.(species)) .* abs.(ks) .* 1e-3) == length(ks)
@test errs[1] < errs[2] < errs[3]
@test errs[1] < tol * 1e-3


# Test an effective low frequency small sphere

    s1 = Specie(
        Acoustic(spatial_dim; ρ=10.2, c=10.1), Sphere(0.001);
        volume_fraction=0.15,
        exclusion_distance = 1.567
    )
    species = [s1];
    basis_field_order = 3

    R = 0.01
    material = Material(Sphere(R),species);

    opts = Dict(
        :tol => tol, :num_wavenumbers => 2,
        :mesh_size => 2.0, :mesh_points => 10,
        :basis_order => basis_order, :basis_field_order => basis_field_order
        , :numberofparticles => material.numberofparticles
    );

    kps = wavenumbers(ωs[1], medium, species; symmetry = PlanarAzimuthalSymmetry(), opts...)

    eff_medium = effective_medium(medium, species; numberofparticles = material.numberofparticles)

    k_low = ωs[1] / eff_medium.c

    @test abs(kps[1] - k_low) / abs(k_low) < 1e-10

    wavemode = WaveMode(ωs[1], kps[1], psource, material;
        basis_order = basis_order,
        basis_field_order = basis_field_order
    )

    # I can not work out theoretically why this scale_number_density should be needed below. However, numerically it clearer is needed.
    scale_number_density = 1.0 - 1.0 / material.numberofparticles
    scat_azi = material_scattering_coefficients(wavemode);

    r = maximum(outer_radius.(species))
    material_lows = [Material(Sphere(R + r1),species) for r1 in (-3.0*r):r/20.0:(r)];

    effective_spheres = map(material_lows) do m
        Particle(eff_medium, m.shape)
    end;

    Linc = basis_field_order + basis_order;

    source_coefficients =  regular_spherical_coefficients(source)(Linc,zeros(3),ωs[1])

    errs =  map(eachindex(material_lows)) do i
        Tmat = MultipleScattering.t_matrix(effective_spheres[i], medium, ωs[1], Linc)
        Tmat = Tmat ./ scale_number_density
        norm(scat_azi - Tmat * source_coefficients) / norm(source_coefficients)
    end

    err, i = findmin(errs)
    @test err < 1e-20
    @test outer_radius(material_lows[i].shape) == R - r

# Test the radius for the scattering coefficients seperate from eigenvectors

    wavemodes = [
        EffectiveRegularWaveMode(ωs[1], kps[1], psource, m, wavemode.eigenvectors;
            basis_order = wavemode.basis_order, basis_field_order = wavemode.basis_field_order)
    for m in material_lows]

    # The scattering coefficients depend on the radius of the material
    scat_azis = material_scattering_coefficients.(wavemodes);

    material_low = Material(Sphere(R - r),species);
    effective_sphere = Particle(eff_medium, material_low.shape);

    Tmat = MultipleScattering.t_matrix(effective_sphere, medium, ωs[1], Linc) ./ scale_number_density

    errs =  map(eachindex(material_lows)) do j
        norm(scat_azis[j] - Tmat * source_coefficients) / norm(source_coefficients)
    end;

    err, j = findmin(errs)
    @test err < 1e-20
    # We expect the min error to be for a materail (which contains all particles) with radius R
    @test outer_radius(material_lows[j].shape) == R

# Test radius for boundary condition seperate from the scattering coefficients

    wavemodes = [
        WaveMode(ωs[1], kps[1], psource, m;
            basis_order = basis_order, basis_field_order = basis_field_order)
    for m in material_lows];

    wavemodes_2 = [
        EffectiveRegularWaveMode(ωs[1], kps[1], psource, material, wavemodes[j].eigenvectors;
            basis_order = wavemodes[j].basis_order, basis_field_order = wavemodes[j].basis_field_order)
    for j in eachindex(material_lows)];

    scat_azis = material_scattering_coefficients.(wavemodes_2);

    errs =  map(eachindex(material_lows)) do j
        norm(scat_azis[j] - Tmat * source_coefficients) / norm(source_coefficients)
    end

    # The boundary conditions for this case are not at all sensitive to small changes in the radius. For this reason, here we do things differently.

    # Let us find the case with the radius R we expect to be accurate
    j = findfirst(R .== [outer_radius(m.shape) for m in material_lows])
    err = errs[j]

    @test err < 1e-25
    @test err < errs[1]
    @test err < errs[end]

end
