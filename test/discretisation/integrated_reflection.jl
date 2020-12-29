using EffectiveWaves, Test
# using OffsetArrays

@testset "Integrated reflection" begin

    # physical parameters
    θin = 0.13
    k=1.; ho = 2

    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    specie = Specie(Particle(
        Acoustic(2; ρ=0.1, c=0.5), ms.Circle(0.6));
        volume_fraction=0.25
    )

    ω = real(k*medium.c)

    # From effective wave theory
    k_eff0 = wavenumber_low_volumefraction(ω, medium, [specie]; basis_order = ho)
    max_x = 15.0*k/imag(k_eff0)
    x = 0.0:0.001:max_x

    normal = [-1.0,0.0] # an outward normal to the surface
    material = Material(Halfspace(normal),[specie])
    source = PlaneSource(medium, [cos(θin),sin(θin)])

    wave0 = WaveMode(ω, k_eff0, source, material; basis_order = ho)

    # wave0 = EffectivePlaneWaveMode(ω, k_eff0, medium, [specie]; θin = θin, basis_order = ho, tol=1e-8)
    wave_avg0 = DiscretePlaneWaveMode(x, wave0, material.shape)

    R = reflection_coefficient(ω, wave_avg0, source, material)
    R_eff = reflection_coefficient(ω, wave0, source, material)

    @test abs(R-R_eff) < 1e-6

    # the matched wave also gives the same reflection coefficient
    m1 = Int(round(length(x)/50))
    m2 = m1 + 100
    wave_avg1 = deepcopy(wave_avg0)
    wave_avg1.amplitudes = wave_avg0.amplitudes[1:m2,:,:]
    wave_avg1.x = wave_avg0.x[1:m2]

    match_wave = MatchPlaneWaveMode{Float64,2}([wave0], wave_avg1, x[m1:m2])
    R_m = reflection_coefficient(ω, match_wave, source, material)

    @test abs(R_m-R_eff)  < 1e-6

    num_wavenumbers = 4
    k_effs = wavenumbers(ω, medium, [specie]; tol = 1e-8,
        basis_order = ho, num_wavenumbers = num_wavenumbers)

    rel_errors = map(k_effs[1:min(5,length(k_effs))]) do k_eff
        wave = WaveMode(ω, k_eff, source, material;basis_order = ho)
        wave_avg = DiscretePlaneWaveMode(x, wave, material.shape)
        R = reflection_coefficient(ω, wave_avg, source, material)
        R_eff = reflection_coefficient(ω, wave, source, material)
        # @test abs(R-R_eff) < 5e-5

        wave_avg1 = deepcopy(wave_avg)
        wave_avg1.amplitudes = wave_avg.amplitudes[1:m2,:,:]
        wave_avg1.x = wave_avg.x[1:m2]

        match_wave = MatchPlaneWaveMode{Float64,2}([wave], wave_avg1, x[m1:m2])
        R_m = reflection_coefficient(ω, match_wave, source, material)
        # @test abs(R_m-R_eff) < 5e-5
        # abs(R_m-R_eff)
        [abs(R - R_eff), abs(R_m - R_eff)]
    end
    @test sum(rel_errors[1] .< 1e-6) == 2
    @test sum(maximum(rel_errors) .< 5e-5) == 2

end
