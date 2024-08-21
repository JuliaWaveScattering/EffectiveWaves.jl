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

    medium = Acoustic(spatial_dim; ρ=1.00001, c=1.00001)
    medium0 = Acoustic(spatial_dim; ρ=1.0, c=1.0)

    # Choose the species
    radius1 = 1.0
    s1 = Specie(
        Acoustic(spatial_dim; ρ=10.0, c=10.0), radius1;
        volume_fraction=0.1
    )
    species = [s1]

    # Define the plate
    normal = [0.0, 0.0, -1.0] # an outward normal to both surfaces of the plate
    Z = 200.0;
    plate = Plate(normal, Z, [0.0, 0.0, Z / 2])

    # define plane wave sources with slightly different material
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

    @test abs(RTs2[1][1] - RTs1[1][1]) / abs(RTs2[1][1]) < 0.02
    @test abs(RTs1[1][2] - RTs2[1][2]) / abs(RTs2[1][2]) < 0.003
    errorsT = [norm([abs(RTs1[i][2] - RTs2[i][2])]) for i in eachindex(ωs)]
    @test norm(errorsT) / norm(norm.(RTs2)) < 0.005
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

@testset "Compare low frequency plate two media" begin

    # Low frequency test
    spatial_dim = 3
    # Choose the frequency
    ωs = 0.0001:0.0002:0.0011
    basis_order = 2

    # Define all properties
    c = 3.0
    ρ = 3.0
    c0 = 1.0
    ρ0 = 1.0
    cs = 10.0
    ρs = 10.0
    ϕ = 0.1

    medium = Acoustic(spatial_dim; ρ = ρ, c = c)
    medium0 = Acoustic(spatial_dim; ρ = ρ0, c = c0)

    # Choose the species
    r = 1.0
    s1 = Specie(
        Acoustic(spatial_dim; ρ = ρs, c = cs), r;
        volume_fraction = ϕ
    );
    species = [s1]

    # Define the plate
    normal = [0.0, 0.0, -1.0] # an outward normal to both surfaces of the plate
    Z = 200.0
    plate = Plate(normal, Z, [0.0, 0.0, Z / 2])

    # define a plane wave source travelling to the material
    source = PlaneSource(medium, [0.0,0.0,1.0])

    # Reflection and transmittion low frequency approximation
    eff_medium = effective_medium(medium0, species)
    amps = [planewave_amplitudes(ω, source, eff_medium, plate) for ω in ωs]
    Ramp = [amps[i][1] for i in eachindex(amps)]
    Tamp = [amps[i][2] for i in eachindex(amps)]
    RTlows = [[Ramp[i], Tamp[i]] for i in eachindex(amps)]

    # Calculate the effective wavenumber and wavemode numerically from the general methods
    kp_arr = [wavenumbers(ω, medium0, species; tol = 1e-6, num_wavenumbers = 2, basis_order = basis_order) for ω in ωs]
    k_effs = [kps[1] for kps in kp_arr]

    material = Material(medium0,plate,species)
    wavemodes = [WaveMode(ωs[i], k_effs[i], source, material; tol = 1e-6, basis_order = basis_order) for i in eachindex(ωs)];

    RTeffs = [reflection_transmission_coefficients(w, source, material) for w in wavemodes]

    @test abs(RTeffs[1][1] - RTlows[1][1]) / abs(RTlows[1][1]) < 0.09
    @test abs(RTeffs[1][2] - RTlows[1][2]) / abs(RTlows[1][2]) < 0.009
    errors = [norm([abs(RTeffs[i][1] - RTlows[i][1]), abs(RTeffs[i][2] - RTlows[i][2])]) for i in eachindex(ωs)]
    @test norm(errors) / norm(norm.(RTlows)) < 0.06
    
end