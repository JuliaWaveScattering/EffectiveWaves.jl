# using OffsetArrays
# include("../src/average_waves/integral_form.jl")
# include("src/average_waves/integral_form.jl")

@testset "Integrated reflection" begin

    # physical parameters
    θin = 0.13
    k=1.; ho = 2
    medium = Medium(1.0,1.0+0.0im)
    ω = real(k*medium.c)
    specie = Specie(ρ=0.1,r=0.1, c=0.5, volfrac=0.1)

    # From effective wave theory
    k_eff0 = wavenumber_low_volfrac(ω, medium, [specie]; tol = 1e-12)
    max_x = 10.*k/imag(k_eff0)
    x = 0.0:0.001:max_x

    wave0 = EffectiveWave(ω, k_eff0, medium, [specie]; θin = θin, hankel_order = ho, tol=1e-8)
    wave_avg0 = AverageWave(x, wave0)

    R = reflection_coefficient(ω, wave_avg0, medium, specie; θin = θin)
    R_eff = reflection_coefficient(ω, wave0, medium, [specie]; θin = θin, hankel_order = ho)

    @test abs(R-R_eff)/abs(R_eff) < 1e-4 #

    k_effs = wavenumbers(ω, medium, [specie]; tol = 1e-8, hankel_order = ho)

    rel_errors = map(k_effs) do k_eff
        wave = EffectiveWave(ω, k_eff, medium, [specie]; θin = θin, hankel_order = ho)
        wave_avg = AverageWave(x, wave)
        R = reflection_coefficient(ω, wave_avg, medium, specie; θin = θin)
        R_eff = reflection_coefficient(ω, wave, medium, [specie]; θin = θin)
        abs(R-R_eff)/abs(R_eff)
    end
    @test rel_errors[1] < 1e-6 && maximum(rel_errors) < 2e-4
end
