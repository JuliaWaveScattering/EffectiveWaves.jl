using EffectiveWaves, Test, LinearAlgebra

@testset "low frequency 3D planar symmetry" begin

   # Low frequency test
    spatial_dim = 3
    # Choose the frequency
    ω = 1e-5

    # medium = Acoustic(spatial_dim; ρ=0.3, c=0.5)
    φ = 0.2
    ρ = 0.3
    ρo = 10.2; co = 5.2 + 6.0im
    Δρ = φ * (ρ - ρo) / (ρ + 2ρo)

    ks = ω * (1 + 20000.0im)
    # below we choose the background wavespeed to obtain the above effective wavenumber ks in the low frequency limit
    c = (ω * co * sqrt(-1 + Δρ + 0.0im) * sqrt(ρo) * sqrt(-1 + φ + 0.0im)) /
        sqrt(co^2 * ks^2 * (1 + 2*Δρ) * ρo + ω^2 * (-1 + Δρ) * ρ * φ  + 0.0im)

    medium = Acoustic(spatial_dim; ρ=ρ, c=c)

    # Choose the species
    radius1 = 0.001
    s1 = Specie(
        Acoustic(spatial_dim; ρ=ρo, c=co), radius1;
        volume_fraction=φ
    );
    species = [s1]

    k = ω / medium.c

    # For the limit of low frequencies we can define
    eff_medium = effective_medium(medium, species)

    @test abs(ks - ω / eff_medium.c) / abs(ks) < 1e-10
    width = 200.0 # plate width
    @test abs(exp(im*ks*width)) < 1e-10

    # Define a plate
    r = maximum(outer_radius.(species))
    normal = [0.0,0.0,-1.0] # an outward normal to both surfaces of the plate
    plate = Plate(normal,width,[0.0,0.0,width/2])

    plate_low = Plate(normal,width - 2r,[0.0,0.0,width/2])
    halfspace = Halfspace(normal)
    halfspace_low = Halfspace(normal,halfspace.origin-r)

    # define a plane wave source travelling at a 45 degree angle in relation to the material
    # source = PlaneSource(medium, [cos(pi/4.0),0.0,sin(pi/4.0)])
    source = PlaneSource(medium, [0.0,0.0,1.0])

    RT1s = reflection_transmission_coefficients(ω, source, eff_medium, halfspace_low)

    amps = planewave_amplitudes(ω, source, eff_medium, plate_low);
    Ramp = amps[1]
    Tamp = amps[2]

    @test abs(Ramp-RT1s[1]) / abs(Ramp) < 1e-10
    @test abs(amps[3]-RT1s[2]) / abs(RT1s[2]) < 1e-10
    @test abs(Tamp) < 1e-10

    # Calculate the effective wavenumber and wavemode numerically from the general methods
    k_effs = wavenumbers(ω, medium, species; tol = 1e-6, num_wavenumbers = 1, basis_order = 1)
    k_eff = k_effs[1]

    @test abs(k_eff - ks) < 1e-7

    material = Material(halfspace,species)

    # Calculate the wavemode for the first wavenumber
    # the WaveMode function calculates the types of waves and solves the needed boundary conditions
    wavemodes = WaveMode(ω, k_eff, source, material; tol = 1e-6, basis_order = 1);

    Reff = reflection_coefficient(wavemodes, source, material)
    # the material shape for the low frqeuency homoegneous method is r smaller than the effective waves method. So a phase correction is needed.
    @test abs(Reff*exp(im*k*r) - Ramp) < 1e-8

    material = Material(plate,species)
    wavemodes = WaveMode(ω, k_eff, source, material; tol = 1e-6, basis_order = 1);
    RTeff = reflection_transmission_coefficients(wavemodes, source, material);

    @test abs(Reff - RTeff[1]) < 1e-15
    @test abs(RTeff[2]*exp(im*k*r) - Tamp) < 1e-15

end

@testset "3D planar compare planar-azimuthal to planar" begin

    spatial_dim = 3
    # Choose the frequency
    ω = 0.3
    φ = 0.2

    medium = Acoustic(spatial_dim; ρ=1.1, c=0.9)

    # Choose the species
    radius1 = 1.0
    s1 = Specie(
        Acoustic(spatial_dim; ρ=6.0, c=10.0), radius1;
        volume_fraction=φ
    );
    species = [s1]

    k = ω / medium.c

    basis_order = 4

    k_effs = wavenumbers(ω, medium, species; tol = 1e-6, num_wavenumbers = 2, basis_order = basis_order)
    k_eff = k_effs[1]

    halfspace = Halfspace([0.0,0.0,-1.0])
    material = Material(halfspace,species)

    eigs = eigenvectors(ω, k_eff, medium, species, PlanarSymmetry{3}(); basis_order = basis_order)
    azi_eigs = eigenvectors(ω, k_eff, medium, species, PlanarAzimuthalSymmetry{3}(); basis_order = basis_order)

    eigs2 = convert_eigenvector_basis(medium,PlanarAzimuthalSymmetry{3}(),azi_eigs)

    @test abs(abs(dot(eigs2,eigs2)) - norm(eigs) * norm(eigs2)) < 1e-15

    psource0 = PlaneSource(medium, [0.0,0.0,1.0])
    psource = PlaneSource(medium, [1e-3,0.0,1.0])
    source = plane_source(medium; direction = [0.0,0.0,1.0])

    Symmetry(psource,material)
    Symmetry(psource0,material)

    pwavemodes = WaveMode(ω, k_eff, psource0, material; tol = 1e-6, basis_order = basis_order);
    wavemodes = WaveMode(ω, k_eff, psource, material; tol = 1e-6, basis_order = basis_order);

    @test maximum(abs.(wavemodes.eigenvectors - pwavemodes.eigenvectors)) < 1e-3

end

@testset "Compare 3D plate and halfspace" begin

   # Low frequency test
    spatial_dim = 3
    # Choose the frequency
    ω = 0.1
    basis_order = 2

    # medium = Acoustic(spatial_dim; ρ=0.3, c=0.5)
    φ = 0.2
    ρ = 0.3
    ρo = 10.2; co = 5.2 + 6.0im
    Δρ = φ * (ρ - ρo) / (ρ + 2ρo)

    ks = ω * (1 + 20.0im)
    # below we choose the background wavespeed to obtain the above effective wavenumber ks in the low frequency limit
    c = (ω * co * sqrt(-1 + Δρ + 0.0im) * sqrt(ρo) * sqrt(-1 + φ + 0.0im)) /
        sqrt(co^2 * ks^2 * (1 + 2*Δρ) * ρo + ω^2 * (-1 + Δρ) * ρ * φ  + 0.0im)

    medium = Acoustic(spatial_dim; ρ=ρ, c=c)

    # Choose the species
    radius1 = 1.0
    s1 = Specie(
        Acoustic(spatial_dim; ρ=ρo, c=co), radius1;
        volume_fraction=φ
    );
    species = [s1]

    k = ω / medium.c

    # For the limit of low frequencies we can define
    eff_medium = effective_medium(medium, species)
    ω / eff_medium.c

    width = 20.0 # plate width
    @test abs(exp(im*ks*width)) < 1e-10

    # Define a plate
    r = maximum(outer_radius.(species))
    normal = [0.0,0.0,-1.0] # an outward normal to both surfaces of the plate
    plate = Plate(normal,width,[0.0,0.0,width/2])
    halfspace = Halfspace(normal)

    # define a plane wave source travelling at a 45 degree angle in relation to the material
    # source = PlaneSource(medium, [cos(pi/4.0),0.0,sin(pi/4.0)])
    source = PlaneSource(medium, [0.0,0.0,1.0])

    # Calculate the effective wavenumber and wavemode numerically from the general methods
    k_effs = wavenumbers(ω, medium, species; tol = 1e-6, num_wavenumbers = 2, basis_order = basis_order)
    k_eff = k_effs[1]

    material = Material(halfspace,species)
    wavemodes = WaveMode(ω, k_eff, source, material; tol = 1e-6, basis_order = basis_order);

    Reff = reflection_coefficient(wavemodes, source, material)

    material = Material(plate,species)
    wavemodes = WaveMode(ω, k_eff, source, material; tol = 1e-6, basis_order = basis_order);
    RTeff = reflection_transmission_coefficients(wavemodes, source, material);

    @test abs(Reff - RTeff[1]) < 1e-10
end

@testset "Compare low frequency halfspace" begin

   # Low frequency test
    spatial_dim = 3
    # Choose the frequency
    ωs = 0.001:0.02:0.1
    basis_order = 2

    medium = Acoustic(spatial_dim; ρ=1.0, c=1.0)

    # Choose the species
    radius1 = 1.0
    s1 = Specie(
        Acoustic(spatial_dim; ρ=10.0, c=10.0), radius1;
        volume_fraction=0.05
    );
    species = [s1]

    # k = ω / medium.c

    # For the limit of low frequencies we can define
    eff_medium = effective_medium(medium, species)
    # ω / eff_medium.c

    # Define the halfspace
    r = maximum(outer_radius.(species))
    normal = [0.0,0.0,-1.0] # an outward normal to both surfaces of the plate
    halfspace = Halfspace(normal)

    # define a plane wave source travelling at a 45 degree angle in relation to the material
    # source = PlaneSource(medium, [cos(pi/4.0),0.0,sin(pi/4.0)])
    source = PlaneSource(medium, [0.0,0.0,1.0])

    # reflection low freq approx
    Rlows = [reflection_transmission_coefficients(ω, source, eff_medium, halfspace)[1] for ω in ωs];

    # Calculate the effective wavenumber and wavemode numerically from the general methods
    kp_arr = [wavenumbers(ω, medium, species; tol = 1e-6, num_wavenumbers = 2, basis_order = basis_order) for ω in ωs]
    k_effs = [kps[1] for kps in kp_arr]

    material = Material(halfspace,species)
    wavemodes = [WaveMode(ωs[i], k_effs[i], source, material; tol = 1e-6, basis_order = basis_order) for i in eachindex(ωs)];

    w1 = WaveMode(ωs[end], k_effs[end], source, material; tol = 1e-6, basis_order = basis_order)
    e1 = w1.eigenvectors / norm(w1.eigenvectors)
    w2 = WaveMode(ωs[end], -k_effs[end], source, material; tol = 1e-6, basis_order = basis_order)
    e2 = w2.eigenvectors / norm(w2.eigenvectors)

    Reffs = [reflection_coefficient(w, source, material) for w in wavemodes]

    @test abs(Reffs[1] - Rlows[1]) / abs(Rlows[1]) < 0.01
    @test norm(Reffs - Rlows) / norm(Rlows) < 0.1
end
