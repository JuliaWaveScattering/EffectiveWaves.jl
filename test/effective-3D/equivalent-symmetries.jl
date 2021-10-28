using EffectiveWaves, Test, LinearAlgebra

@testset "symmetry equivalent wavenumbers and eigenvectors" begin

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

s1 = Specie(
    Acoustic(spatial_dim; ρ=0.2, c=0.1), Sphere(0.4);
    volume_fraction=0.2
);
s2 = Specie(
    Acoustic(spatial_dim; ρ=10.2, c=10.1), Sphere(0.2);
    volume_fraction=0.1
);
species = [s1,s2]
# species = [s1]

ω = 0.9
tol = 1e-7

basis_order = 2

### Test the equivalence between dispersion equations

    R_det = dispersion_equation(ω, medium, species, RadialSymmetry{spatial_dim}(); basis_order = basis_order)
    AP_det = dispersion_equation(ω, medium, species, PlanarAzimuthalSymmetry{spatial_dim}(); basis_order = basis_order)
    P_det = dispersion_equation(ω, medium, species, PlanarSymmetry{spatial_dim}(); basis_order = basis_order)
    AR_det = dispersion_equation(ω, medium, species, AzimuthalSymmetry{spatial_dim}(); basis_order = basis_order)
    Reg_det = dispersion_equation(ω, medium, species, WithoutSymmetry{spatial_dim}(); basis_order = basis_order)

    # Lighter to calculate wavenumbers (also known as the eigenvalues) when assume more symmmetries

    # The RadialSymmetry eigensystem is almost identical to the PlanarAzimuthalSymmetry system
    # R_kps = wavenumbers(ω, medium, species;
    #     num_wavenumbers = 2, tol = tol,  basis_order = basis_order,
    #     symmetry = RadialSymmetry{3}())

    AP_kps = wavenumbers(ω, medium, species;
        num_wavenumbers = 2, tol = tol,  basis_order = basis_order,
        symmetry = PlanarAzimuthalSymmetry())

    P_kps = wavenumbers(ω, medium, species;
        num_wavenumbers = 2, tol = tol, basis_order = basis_order,
        symmetry = PlanarSymmetry{3}())

    # The dispersion equations R_det become very unstable when using effective wavenumbers with higher dispersion.
    AP_kps = AP_kps[1:min(length(AP_kps),6)]
    P_kps = P_kps[1:min(length(P_kps),6)]

    is = findall(AP_det.(P_kps) .< tol)
    AP_kps = sort(reduce_kvecs([AP_kps; P_kps[is]],tol), by=imag)

    # As plane waves with azimuthal symmetry is a sub-case of plane-waves, and all materials allow for the effective wavenumbers of plane waves, all the below determinant equations should be satisfied

    @test maximum(AP_det.(AP_kps)) < tol
    # @test maximum(AP_det.(R_kps)) < tol
    @test maximum(R_det.(AP_kps)) < tol
    @test maximum(P_det.(AP_kps)) < tol
    @test maximum(AR_det.(AP_kps)) < tol^2
    @test maximum(Reg_det.(AP_kps)) < tol^3

    P_disp = dispersion_complex(ω, medium, species,  PlanarSymmetry{spatial_dim}(); tol = tol, basis_order = basis_order)

    # However, there do exist effective wavenumbers for plane-waves which have eigen-vectors that do not satisfy azimuthal symmetry. This is why maximum(AP_det.(P_kps)) != 0.0
    @test maximum(P_det.(P_kps)) < tol
    @test maximum(AR_det.(P_kps)) < tol^2
    @test maximum(Reg_det.(P_kps)) < tol^3

### Test that any regular eigenvector can be transform into an plane-wave eigenvector

    basis_order = 2
    basis_field_order = 2*basis_order

    Reg_MM = eigensystem(ω, medium, species, WithoutSymmetry{spatial_dim}();
            basis_order = basis_order,
            basis_field_order = basis_field_order
    )

    k_eff = P_kps[1]
    MM_svd = svd(Reg_MM(k_eff))
    inds = findall(abs.(MM_svd.S) .< 1e-4)

    # Rvs = [MM_svd.V[:,i] for i in inds] # eigenvectors
    # RvMs = [transpose(reshape(v, :, (basis_order+1)^2)) for v in Rvs]
    # Pvs = [
    #     sum([ vM[i] * Ys[i[2]] / 1.0im^ls[i[2]] for i in CartesianIndices(vM)], dims=2)[:]
    # for vM in RvMs]


    RvM = MM_svd.V[:,inds]
    S = length(species)
    RvM = reshape(RvM,(:,S,size(RvM,2)))

    # θp = 0.2; φp = 0.1;
    direction_eff = [0.2,-0.5, 1.0]
    rθφ_eff = cartesian_to_radial_coordinates(direction_eff)

    Ys = spherical_harmonics(basis_field_order, rθφ_eff[2], rθφ_eff[3]);
    ls, ms = spherical_harmonics_indices(basis_field_order)

    L = basis_order
    L1 = basis_field_order

    RvM = reshape(RvM,(:,(basis_order+1)^2,S,size(RvM)[end]))
    Pvs = reshape(sum([RvM[i] * Ys[i[1]] / 1.0im^ls[i[1]] for i in CartesianIndices(RvM)], dims=1),(:,size(RvM)[end]))

    P_MM = eigensystem(ω, medium, species, PlanarSymmetry{spatial_dim}();
            direction_eff = direction_eff,
            basis_order = basis_order
    )

    PM = P_MM(k_eff)

    @test maximum(norm(PM * Pvs[:,i]) for i in axes(Pvs,2)) < 1e-8

end

@testset "Equivalence for solving sphere scattering" begin

### Test that the complete regular solutions is the same as the azimuthal solution when there is azimuthal symmetry

    spatial_dim = 3
    medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

    s1 = Specie(
        Acoustic(spatial_dim; ρ=0.1, c=0.1), Sphere(0.4);
        volume_fraction=0.2
    );
    s2 = Specie(
        Acoustic(spatial_dim; ρ=7.2, c=10.1), Sphere(0.3);
        volume_fraction=0.1
    );
    species = [s1,s2]
    # species = [s1]

    ω = 0.9
    tol = 1e-7

    basis_order = 2
    basis_field_order = 2*basis_order

    kps = wavenumbers(ω, medium, species;
        num_wavenumbers = 2, tol = tol,  basis_order = basis_order,
        symmetry = PlanarAzimuthalSymmetry())
    k_eff = kps[1]

    θ = 0.0
    psource = PlaneSource(medium, [sin(θ),0.0,cos(θ)]);

    # change direction slightly to break symmetry
    source = plane_source(medium; direction = [sin(θ),1e-7,cos(θ)])
    Symmetry(source) == PlanarSymmetry{3}()

    material = Material(Sphere(4.0),species);

    nn1_indexs = [
            [l,m,l1,m1]
        for l = 0:basis_order for m = -l:l
    for l1 = 0:basis_field_order for m1 = -l1:l1];

    # test boundary conditions for all wavenumbers
    # [
    #     eigenvectors(ω, kp, source, material;
    #         basis_order = basis_order,
    #         basis_field_order = basis_field_order)
    # for kp in AP_kps];
    #
    # [
    #     eigenvectors(ω, kp, psource, material;
    #         basis_order = basis_order,
    #         basis_field_order = basis_field_order)
    # for kp in AP_kps];

    wave_azi = WaveMode(ω, k_eff, psource, material;
        basis_order = basis_order, basis_field_order = basis_field_order)

    wave_reg = WaveMode(ω, k_eff, source, material;
        basis_order = basis_order, basis_field_order = basis_field_order)

    A0inds = findall([
            m1 != -m
        for l = 0:basis_order for m = -l:l
    for l1 = 0:basis_field_order for m1 = -l1:l1]);

    @test maximum(abs.(sum(wave_reg.eigenvectors,dims=(2,3))[A0inds])) < tol

    Ainds = findall([
            m1 == -m
        for l = 0:basis_order for m = -l:l
    for l1 = 0:basis_field_order for m1 = -l1:l1]);

    azi_vecs = sum(wave_azi.eigenvectors, dims=(2,3))[:]
    reg_azi_vecs = sum(wave_reg.eigenvectors, dims=(2,3))[Ainds]

    @test norm(reg_azi_vecs - azi_vecs) / norm(azi_vecs) < tol

    ### Check that the predicting average scattering from a sphere and an incident plane-wave is the same

    scat_azi = material_scattering_coefficients(wave_azi);
    scat_reg = material_scattering_coefficients(wave_reg);

    @test norm(scat_azi - scat_reg) / norm(scat_azi) < tol

### Test that the radially symmetric solution is the same as the complete regular solutions and the azimuthal solution when there is radial symmetry

    source = regular_spherical_source(medium, [1.0+0.0im];
       position = [0.0,0.0,0.0], symmetry = RadialSymmetry{3}()
    );

    # Same source but pretend it has no symmetry! These tests are getting weirder and weirder..
    source_reg = regular_spherical_source(medium, [1.0+0.0im];
       position = [0.0,0.0,0.0], symmetry = WithoutSymmetry{3}()
    );

    @test Symmetry(source,material) == RadialSymmetry{3}()

    eigs_rad = eigenvectors(ω, k_eff, source, material; basis_order = basis_order)
    α = solve_boundary_condition(ω, k_eff, eigs_rad, source, material;
            basis_order = basis_order)

    wave_reg = WaveMode(ω, k_eff, source_reg, material;
        basis_order = basis_order, basis_field_order = basis_field_order)

    eigs_reg = sum(wave_reg.eigenvectors,dims=(3))

    rad_inds = findall([
            m1 == 0 && m == 0 && l == l1
        for l = 0:basis_order for m = -l:l
    for l1 = 0:basis_field_order for m1 = -l1:l1]);

    # rad_inds = findall([
    #         m1 == - m && l == l1
    #     for l = 0:basis_order for m = -l:l
    # for l1 = 0:basis_field_order for m1 = -l1:l1])

    findall(abs.(eigs_reg) .> 1e-12)

    eigs_reg[rad_inds,:]

    eigs_rad .*  (eigs_reg[1,1] / eigs_rad[1,1])

    eigs_reg[rad_inds,:]
    eigs_rad .* α[1]

    sum(eigs_reg[rad_inds,:],dims=2)
    sum(eigs_rad .* α[1],dims=2)

    # Calculate the eigenvectors
    radwave = WaveMode(ω, k_eff, source, material;
        basis_order = basis_order,
    )

end
# inds = findall(waves .!= nothing)
#
# ωs = ωs[inds]
# kas = kas[inds]
# kps = kps[inds]
# waves = waves[inds]
# basis_orders = basis_orders[inds]
# basis_field_orders = basis_field_orders[inds]
#
# scat_coefs = material_scattering_coefficients.(waves);
# # coef = material_scattering_coefficients(waves[i]);
# # coef3 = material_scattering_coefficients(w);
#
# # Calculate low frequency scattering
#
# material_low = Material(
#     Sphere(outer_radius(material.shape) - r),
#     species
# );
#
# effective_sphere = Particle(eff_medium, material_low.shape);
#
# scat_coef_lows = map(eachindex(ωs)) do i
#     Linc = basis_field_orders[i] + basis_orders[i]
#     source_coefficients = regular_spherical_coefficients(source)(Linc,zeros(3),ωs[i]);
#     Tmat = t_matrix(effective_sphere, medium, ωs[i], Linc);
#     Tmat * source_coefficients
# end
