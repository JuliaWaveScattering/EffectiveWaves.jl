using EffectiveWaves, Test
# using OffsetArrays

@testset "Integrated reflection" begin

    # physical parameters
    θin = 0.13
    k=1.; ho = 3
    medium = Medium(1.0,1.0+0.0im)
    ω = real(k*medium.c)
    specie = Specie(ρ=0.1,r=0.6, c=0.5, volfrac=0.1)

    # From effective wave theory
    k_eff0 = wavenumber_low_volfrac(ω, medium, [specie]; tol = 1e-12)
    max_x = 12.0*k/imag(k_eff0)
    x = 0.0:0.002:max_x

    wave0 = EffectiveWave(ω, k_eff0, medium, [specie]; θin = θin, basis_order = ho, tol=1e-8)
    wave_avg0 = AverageWave(x, wave0)

    R = reflection_coefficient(ω, wave_avg0, medium, specie; θin = θin)
    R_eff = reflection_coefficient(ω, wave0, medium, [specie]; θin = θin, basis_order = ho)

    @test abs(R-R_eff) < 2e-6 #

    # the matched wave also gives the same reflection coefficient
    m1 = Int(round(length(x)/50))
    m2 = m1 + 100
    wave_avg1 = deepcopy(wave_avg0)
    wave_avg1.amplitudes = wave_avg0.amplitudes[1:m2,:,:]
    wave_avg1.x = wave_avg0.x[1:m2]

    match_wave = MatchWave{Float64}([wave0], wave_avg1, x[m1:m2])
    R_m = reflection_coefficient(ω, match_wave, medium, specie; θin = θin)
    @test abs(R_m-R_eff)  < 1e-6

    num_wavenumbers = 4
    k_effs = wavenumbers(ω, medium, [specie]; tol = 1e-8,
        basis_order = ho, num_wavenumbers = num_wavenumbers)

    rel_errors = map(k_effs[1:end]) do k_eff
        wave = EffectiveWave(ω, k_eff, medium, [specie]; θin = θin, basis_order = ho)
        wave_avg = AverageWave(x, wave)
        R = reflection_coefficient(ω, wave_avg, medium, specie; θin = θin)
        R_eff = reflection_coefficient(ω, wave, medium, [specie]; θin = θin)
        @test abs(R-R_eff) < 5e-5

        wave_avg1 = deepcopy(wave_avg)
        wave_avg1.amplitudes = wave_avg.amplitudes[1:m2,:,:]
        wave_avg1.x = wave_avg.x[1:m2]

        match_wave = MatchWave{Float64}([wave], wave_avg1, x[m1:m2])
        R_m = reflection_coefficient(ω, match_wave, medium, specie; θin = θin)
        @test abs(R_m-R_eff) < 5e-5
    end

end
