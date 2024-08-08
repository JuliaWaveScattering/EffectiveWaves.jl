using EffectiveWaves, Test, LinearAlgebra

@testset "Compare one and two media halfspace" begin

    # Two media test
    spatial_dim = 3
    # Choose the frequency
    ωs = 0.001:0.5:3.001
    basis_order = 2

    medium = Acoustic(spatial_dim; ρ=1.001, c=1.001)
    medium0 = Acoustic(spatial_dim; ρ=1.0, c=1.0)

    # Choose the species
    radius1 = 1.0
    s1 = Specie(
        Acoustic(spatial_dim; ρ=10.0, c=10.0), radius1;
        volume_fraction=0.2
    )
    species = [s1]

    # Define the halfspace
    r = maximum(outer_radius.(species))
    normal = [0.0, 0.0, -1.0] # an outward normal to both surfaces of the plate
    halfspace = Halfspace(normal)

    # define a plane wave sources with slightly different material
    source1 = PlaneSource(medium0, [0.0, 0.0, 1.0])
    source2 = PlaneSource(medium, [0.0, 0.0, 1.0])

    # Calculate the effective wavenumber and wavemode numerically from the general methods
    kp_arr = [wavenumbers(ω, medium0, species; tol=1e-6, num_wavenumbers=2, basis_order=basis_order) for ω in ωs]
    k_effs = [kps[1] for kps in kp_arr]

    # Define material 
    material = Material(medium0, halfspace, species)

    # Compute both wavemodes, for one and two slightly different media
    wavemodes1 = [WaveMode(ωs[i], k_effs[i], source1, material; tol=1e-6, basis_order=basis_order) for i in eachindex(ωs)]
    wavemodes2 = [WaveMode(ωs[i], k_effs[i], source2, material; tol=1e-6, basis_order=basis_order) for i in eachindex(ωs)]

    # Compute reflection coefficients for both cases
    Reffs1 = [reflection_coefficient(w, source1, material) for w in wavemodes1]
    Reffs2 = [reflection_coefficient(w, source2, material) for w in wavemodes2]

    @test abs(Reffs2[1] - Reffs1[1]) / abs(Reffs1[1]) < 0.01
    @test norm(Reffs2 - Reffs1) / norm(Reffs1) < 0.01
end

@testset "Compare one and two media plate" begin

    # Two media test
    spatial_dim = 3
    # Choose the frequency
    ωs = 0.001:0.5:3.001
    basis_order = 2

    medium = Acoustic(spatial_dim; ρ=1.001, c=1.001)
    medium0 = Acoustic(spatial_dim; ρ=1.0, c=1.0)

    # Choose the species
    radius1 = 1.0
    s1 = Specie(
        Acoustic(spatial_dim; ρ=0.01, c=0.01), radius1;
        volume_fraction=0.2
    )
    species = [s1]

    # Define the halfspace
    r = maximum(outer_radius.(species))
    normal = [0.0, 0.0, -1.0] # an outward normal to both surfaces of the plate
    width = 200.0;
    plate = Plate(normal, width, [0.0, 0.0, width / 2])

    # define a plane wave sources with slightly different material
    source1 = PlaneSource(medium0, [0.0, 0.0, 1.0])
    source2 = PlaneSource(medium, [0.0, 0.0, 1.0])

    # Calculate the effective wavenumber and wavemode numerically from the general methods
    kp_arr = [wavenumbers(ω, medium0, species; tol=1e-6, num_wavenumbers=2, basis_order=basis_order) for ω in ωs]
    k_effs = [kps[1] for kps in kp_arr]

    # Define material 
    material = Material(medium0, plate, species)

    # Compute both wavemodes, for one and two slightly different media
    wavemodes1 = [WaveMode(ωs[i], k_effs[i], source1, material; tol=1e-6, basis_order=basis_order) for i in eachindex(ωs)]
    wavemodes2 = [WaveMode(ωs[i], k_effs[i], source2, material; tol=1e-6, basis_order=basis_order) for i in eachindex(ωs)]

    # Compute reflection coefficients for both cases
    RTs1 = [reflection_transmission_coefficients(w, source1, material) for w in wavemodes1]
    RTs2 = [reflection_transmission_coefficients(w, source2, material) for w in wavemodes2]

    @test norm(RTs2[1] - RTs1[1]) / norm(RTs1[1]) < 0.001
    @test sum(norm.(RTs2 - RTs1)) / sum(norm.(RTs1)) < 0.01
end

@testset "Compare low frequency halfspace two media" begin

    # Low frequency test
    spatial_dim = 3
    # Choose the frequency
    ωs = 0.0001:0.0002:0.0011
    basis_order = 2

    medium = Acoustic(spatial_dim; ρ = 3.0, c = 3.0)
    medium0 = Acoustic(spatial_dim; ρ = 1.0, c = 1.0)

    # Choose the species
    radius1 = 1.0
    s1 = Specie(
        Acoustic(spatial_dim; ρ=10.0, c=10.0), radius1;
        volume_fraction=0.2
    );
    species = [s1]

    # For the limit of low frequencies we can define
    eff_medium = effective_medium(medium0, species)

    # Define the halfspace
    r = maximum(outer_radius.(species))
    normal = [0.0,0.0,-1.0] # an outward normal to both surfaces of the plate
    halfspace = Halfspace(normal)

    # define a plane wave source travelling at a 45 degree angle in relation to the material
    # source = PlaneSource(medium, [cos(pi/4.0),0.0,sin(pi/4.0)])
    source = PlaneSource(medium, [0.0,0.0,1.0])

    # reflection low freq approx
    Rlows = [reflection_coefficient(ω, source, eff_medium, halfspace) for ω in ωs];

    # Calculate the effective wavenumber and wavemode numerically from the general methods
    kp_arr = [wavenumbers(ω, medium0, species; tol = 1e-6, num_wavenumbers = 2, basis_order = basis_order) for ω in ωs]
    k_effs = [kps[1] for kps in kp_arr]

    material = Material(medium0,halfspace,species)
    wavemodes = [WaveMode(ωs[i], k_effs[i], source, material; tol = 1e-6, basis_order = basis_order) for i in eachindex(ωs)];

    Reffs = [reflection_coefficient(w, source, material) for w in wavemodes]

    @test abs(Reffs[1] - Rlows[1]) / abs(Rlows[1]) < 0.00002
    @test norm(Reffs - Rlows) / norm(Rlows) < 0.0002
end


