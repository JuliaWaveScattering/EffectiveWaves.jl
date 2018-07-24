# using OffsetArrays
# include("../src/average_waves/integral_form.jl")
# include("src/average_waves/integral_form.jl")

@testset "Integrated reflection" begin

    # physical parameters
    θin = 0.0
    k=1.; ho = 2
    medium = Medium(1.0,1.0+0.0im)
    ω = real(k*medium.c)
    specie = Specie(ρ=0.1,r=0.1, c=0.5, volfrac=0.1)
    specie2 = Specie(ρ=2.0, r=1.0, c=0.1, volfrac=0.15)

    # From effective wave theory
    k_eff0 = wavenumber_low_volfrac(ω, medium, [specie]; tol = 1e-12)
    max_x = 10.*k/imag(k_eff0)
    x = 0.0:0.001:max_x

    amps0_eff = AverageWave(ω, x, medium, [specie];
            k_eff = k_eff0, tol = 1e-12, θin=θin, hankel_order = ho)

    R = reflection_coefficient_integrated(ω, medium, specie; amps = amps0_eff, θin = θin, hankel_order = ho)

    wave0 = EffectiveWave(ω, k_eff0, medium, [specie]; θin = θin, hankel_order = ho, tol=1e-8)
    R_eff = reflection_coefficient(ω, wave0, medium, [specie]; θin = θin, hankel_order = ho)

    @test abs(R-R_eff)/abs(R_eff) < 1e-4 #

    k_effs = wavenumbers(ω, medium, [specie]; mesh_points = 10, tol = 1e-8, hankel_order = ho)
    k_effs = sort(k_effs, by=imag)

    k_effs2 = wavenumbers(ω, medium, [specie2]; tol = 1e-8, hankel_order = ho)
    k_effs2 = sort(k_effs2, by=imag)

    amps_eff1 = AverageWave(ω, x, medium, [specie];
            k_eff = k_effs[1], hankel_order = ho, θin=θin, tol=1e-8)
    amps_eff2 = AverageWave(ω, x, medium, [specie];
            k_eff = k_effs[2], hankel_order = ho, θin=θin, tol=1e-8)
    amps_eff3 = AverageWave(ω, x, medium, [specie];
            k_eff = k_effs[end], hankel_order = ho, θin=θin, tol=1e-8)
    amps2_eff = AverageWave(ω, x, medium, [specie2];
            k_eff = k_effs2[1], hankel_order = ho, θin=θin, tol=1e-8)

    R = reflection_coefficient_integrated(ω, medium, specie; amps = amps_eff1, θin = θin, hankel_order = ho)
    wave = EffectiveWave(ω, k_effs[1], medium, [specie]; θin = θin, hankel_order = ho)
    R_eff = reflection_coefficient(ω, wave, medium, [specie]; θin = θin)
    @test abs(R-R_eff)/abs(R_eff) < 1e-6

    R = reflection_coefficient_integrated(ω, medium, specie; amps = amps_eff2, θin = θin, hankel_order = ho)
    wave = EffectiveWave(ω, k_effs[2], medium, [specie]; θin = θin, hankel_order = ho)
    R_eff = reflection_coefficient(ω, wave, medium, [specie]; θin = θin)
    @test abs(R-R_eff)/abs(R_eff) < 5e-5 # the larger the eff wavenumbers, the bigger the integral error

    R = reflection_coefficient_integrated(ω, medium, specie; amps = amps_eff3, θin = θin, hankel_order = ho)
    wave = EffectiveWave(ω, k_effs[end], medium, [specie]; θin = θin, hankel_order = ho)
    R_eff = reflection_coefficient(ω, wave, medium, [specie]; θin = θin)
    @test abs(R-R_eff)/abs(R_eff) < 2e-4

    R = reflection_coefficient_integrated(ω, medium, specie2; amps = amps2_eff, θin = θin, hankel_order = ho)
    wave = EffectiveWave(ω, k_effs2[1], medium, [specie2]; θin = θin, hankel_order = ho)
    R_eff = reflection_coefficient(ω, wave, medium, [specie2]; θin = θin)
    @test abs(R-R_eff)/abs(R_eff) < 1e-5
end
