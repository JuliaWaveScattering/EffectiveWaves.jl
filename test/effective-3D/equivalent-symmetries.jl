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

### Test the equivalence between dispersion equations

    AP_det = dispersion_equation(ω, medium, species, PlanarAzimuthalSymmetry{spatial_dim}())
    P_det = dispersion_equation(ω, medium, species, PlanarSymmetry{spatial_dim}())
    AR_det = dispersion_equation(ω, medium, species, AzimuthalSymmetry{spatial_dim}())
    R_det = dispersion_equation(ω, medium, species, WithoutSymmetry{spatial_dim}())

    # Lighter to calculate wavenumbers (also known as the eigenvalues) when assume more symmmetries
    AP_kps = wavenumbers(ω, medium, species;
        num_wavenumbers = 2, tol = tol,
        symmetry = PlanarAzimuthalSymmetry())

    P_kps = wavenumbers(ω, medium, species;
        num_wavenumbers = 2, tol = tol,
        symmetry = PlanarSymmetry{3}())

    # The dispersion equations R_det become very unstable when using effective wavenumbers with higher dispersion.
    AP_kps = AP_kps[1:min(length(AP_kps),6)]
    P_kps = P_kps[1:min(length(P_kps),6)]

    is = findall(AP_det.(P_kps) .< tol)
    AP_kps = sort(reduce_kvecs([AP_kps; P_kps[is]],tol), by=imag)

    # As plane waves with azimuthal symmetry is a sub-case of plane-waves, and all materials allow for the effective wavenumbers of plane waves, all the below determinant equations should be satisfied
    @test maximum(AP_det.(AP_kps)) < tol
    @test maximum(P_det.(AP_kps)) < tol
    @test maximum(AR_det.(AP_kps)) < tol^2
    @test maximum(R_det.(AP_kps)) < tol^3

    P_disp = dispersion_complex(ω, medium, species,  PlanarSymmetry{spatial_dim}(); tol = tol)

    # However, there do exist effective wavenumbers for plane-waves which have eigen-vectors that do not satisfy azimuthal symmetry. This is why maximum(AP_det.(P_kps)) != 0.0
    @test maximum(P_det.(P_kps)) < tol
    @test maximum(AR_det.(P_kps)) < tol^2
    @test maximum(R_det.(P_kps)) < tol^3

### Test that any regular eigenvector can be transform into an plane-wave eigenvector

    basis_order = 2
    basis_field_order = 2*basis_order

    R_MM = eigensystem(ω, medium, species, WithoutSymmetry{spatial_dim}();
            basis_order = basis_order,
            basis_field_order = basis_field_order
    )

    k_eff = P_kps[1]
    MM_svd = svd(R_MM(k_eff))
    inds = findall(abs.(MM_svd.S) .< 1e-4)

    # Rvs = [MM_svd.V[:,i] for i in inds] # eigenvectors
    # RvMs = [transpose(reshape(v, :, (basis_order+1)^2)) for v in Rvs]
    # Pvs = [
    #     sum([ vM[i] * Ys[i[2]] / 1.0im^ls[i[2]] for i in CartesianIndices(vM)], dims=2)[:]
    # for vM in RvMs]


    RvM = MM_svd.V[:,inds]
    S = length(species)
    RvM = reshape(RvM,(:,S,size(RvM,2)))

    θp = 0.2; φp = 0.1;

    Ys = spherical_harmonics(basis_field_order, θp, φp);
    ls, ms = spherical_harmonics_indices(basis_field_order)

    L = basis_order
    L1 = basis_field_order

    RvM = reshape(RvM,(:,(basis_order+1)^2,S,size(RvM)[end]))
    Pvs = reshape(sum([RvM[i] * Ys[i[1]] / 1.0im^ls[i[1]] for i in CartesianIndices(RvM)], dims=1),(:,size(RvM)[end]))


    P_MM = eigensystem(ω, medium, species, PlanarSymmetry{spatial_dim}();
            θp = θp, φp = φp,
            basis_order = basis_order
    )

    PM = P_MM(k_eff)

    @test maximum(norm(PM * Pvs[:,i]) for i in axes(Pvs,2)) < 1e-8

### Test that the complete regular solutions is the same as the azimuthal solution hen there is azimuthal symmetry

    θ = 0.0
    psource = PlaneSource(medium, [sin(θ),0.0,cos(θ)]);
    source = plane_source(medium; direction = [sin(θ),0.0,cos(θ)])
    material = Material(Sphere(4.0),species);

    nn1_indexs = [
            [l,m,l1,m1]
        for l = 0:basis_order for m = -l:l
    for l1 = 0:basis_field_order for m1 = -l1:l1];

    # test boundary conditions for all wavenumbers
    [
        eigenvectors(ω, kp, source, material;
            basis_order = basis_order,
            basis_field_order = basis_field_order)
    for kp in AP_kps];

    [
        eigenvectors(ω, kp, psource, material;
            basis_order = basis_order,
            basis_field_order = basis_field_order)
    for kp in AP_kps];

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

end
