import Base.Test: @testset, @test, @test_throws

using EffectiveWaves

@testset "Summary" begin

@testset "Examples from: Reflection from multi-species.., Proc.R.Soc.(2018)" begin

    include("../examples/concrete/concrete_species.jl")
    include("../examples/concrete/concrete_species_volfrac.jl")
    # include("../examples/concrete/concrete_species_large-freq.jl") # takes longer

    # Takes too long
    include("../examples/emulsion/fluid_species.jl")
    include("../examples/emulsion/fluid_species_volfrac.jl")
    # include("../examples/emulsion/fluid_species_large-freq.jl") # takes longer

    @test true
end

@testset "Tests for the integral form" begin
    include("integral_form_tests.jl")
    test_integrand()

end

@testset "Particle types and constructors" begin
    ω=0.4
    medium = Medium(ρ=10.,c=2.0+3.0im)
    p_dirichlet = Specie(0.0,0.1) # ρ=0.0, r =0.1
    p_neumann = Specie(Inf,0.1) # ρ=0.0, r =0.1
    Z_dirichlet = Zn(ω, p_dirichlet, medium, 0)
    Z_neumann   = Zn(ω, p_neumann, medium, 0)
    @test Z_dirichlet == besselj(0, 0.1*ω/medium.c)/hankelh1(0, 0.1*ω/medium.c)
    @test Z_neumann == besselj(1, 0.1*ω/medium.c)/hankelh1(1, 0.1*ω/medium.c)

    medium = Medium(ρ=0.,c=2.0+3.0im)
    @test_throws ErrorException Zn(ω, p_dirichlet, medium, 0)
    @test Z_neumann == Zn(ω, p_neumann, medium, 0)

    f_tmp(x,y; kws...) = 0.0
    @test gray_square(rand(10),rand(10), f_tmp) == 0.0
end

@testset "weak scatterers" begin
    # background medium
    medium = Medium(1.0,1.0+0.0im)

    # angular frequencies
    ωs = 0.001:5.0:21.

    species = [
        Specie(ρ=10.,r=0.01, c=12., volfrac=0.05),
        Specie(ρ=3., r=0.2, c=2.0, volfrac=0.04)
    ]

    # wavenumbers
    eff_medium = effective_medium(medium, species)

    k_eff_lows = ωs./eff_medium.c
    k_eff_φs = wavenumber_low_volfrac(ωs, medium, species)

    k_effs = wavenumber(ωs, medium, species)

    @test norm(k_effs - k_eff_φs)/norm(k_effs) < 0.0002
    @test norm(k_effs - k_eff_lows)/norm(k_effs) < 0.01
    @test norm(k_effs[1] - k_eff_lows[1])/norm(k_effs[1]) < 1e-7

    # reflection coefficient
    Rs2 = reflection_coefficient(ωs, medium, species)
    # if the wavenumber k_effs are already calculated, the below is faster
    Rs = reflection_coefficient(ωs, k_effs, medium, species)

    @test norm(Rs - Rs2)/norm(Rs) ≈ 0.0

    # Direct incidence
    R_low = reflection_coefficient_halfspace(medium, eff_medium)
    Rs_φs = reflection_coefficient(ωs, k_eff_φs, medium, species)
    # the below takes a low-volfrac expansion for both the wavenumber and reflection coefficient
    Rs_φs2 = reflection_coefficient_low_volfrac(ωs, medium, species)

    len = length(ωs)
    @test norm(Rs_φs - Rs)/len < 1e-5 # already relative to incident wave amplitude = 1
    @test norm(Rs_φs2 - Rs)/len < 1e-4
    @test abs(R_low - Rs[1]) < 1e-7

    # Vary angle of incidence θin
    θs = 0.1:0.3:(π/2)
    R_low = [reflection_coefficient_halfspace(medium, eff_medium; θin = θ) for θ in θs]
    Rs = [reflection_coefficient(ωs, k_effs, medium, species; θin = θ, hankel_order =7) for θ in θs];
    Rs_φs = [reflection_coefficient_low_volfrac(ωs, medium, species; θin = θ, hankel_order =7) for θ in θs];

    @test maximum(abs(R_low[i] - Rs[i][1]) for i in 1:length(R_low)) < 1e-7
    @test maximum(norm(R)/len for R in (Rs_φs - Rs)) < 0.01

end

@testset "high frequency" begin
    medium = Medium(1.0,1.0+0.0im)
    # Large weak scatterers with low volume fraciton
    species = [
        Specie(ρ=10.,r=1.9, c=12., volfrac=0.04),
        Specie(ρ=3., r=0.7, c=2.0, volfrac=0.02)
    ]
    ωs2 = 20.:30.:121

    k_eff_φs = wavenumber_low_volfrac(ωs2, medium, species; tol=1e-5)
    k_effs = wavenumber(ωs2, medium, species; tol=1e-7) # lowering tol to speed up calculation

    @test norm(k_effs - k_eff_φs)/norm(k_effs) < 1e-4

    Rs = reflection_coefficient(ωs2, k_effs, medium, species)
    Rs_φs = reflection_coefficient(ωs2, k_eff_φs, medium, species)
    Rs_φs2 = reflection_coefficient_low_volfrac(ωs2, medium, species)

    len = length(ωs2)
    @test norm(Rs_φs - Rs)/len < 2e-7
    @test norm(Rs_φs2 - Rs)/len < 1e-5 # the incident wave has amplitude 1, so this is a tiny difference
end

# large volume fraction scatterers,  small size amd on strong scatterer. This is a problamatic case.
@testset "large volume fraction and low frequency" begin
    medium = Medium(1.0,1.0+0.0im)
    species = [
        Specie(ρ=5.,r=0.004, c=1.2, volfrac=0.4),
        Specie(ρ=0.3, r=0.002, c=0.4, volfrac=0.3)
    ]
    # angular frequencies
    ωs = 0.001:5.0:21.

    eff_medium = effective_medium(medium, species)
    k_eff_lows = ωs./eff_medium.c
    k_effs = wavenumber(ωs, medium, species)

    @test norm(k_effs - k_eff_lows)/norm(k_effs) < 0.01
    @test norm(k_effs[1] - k_eff_lows[1])/norm(k_eff_lows[1]) < 1e-7

    Rs = reflection_coefficient(ωs, k_effs, medium, species)
    R_low = reflection_coefficient_halfspace(medium, eff_medium)
    R_low2 = reflection_coefficient(ωs, k_eff_lows, medium, species)

    @test norm(R_low - Rs[1]) < 1e-7
    @test norm(R_low2[1] - Rs[1]) < 1e-7
    @test norm(R_low2 - Rs) < 0.005
end

# This case is numerically challenging, because wavenumber() has many roots close together. Make sure spacing in ωs is small to help the optimisation method
@testset "strong scatterers and low frequency" begin
    medium = Medium(1.0,1.0+0.0im)
    species = [
        Specie(ρ=5.,r=0.004, c=0.002, volfrac=0.2),
        Specie(ρ=0.3, r=0.002, c=0.01, volfrac=0.1)
    ]

    ωs = 0.001:0.002:0.01
    eff_medium = effective_medium(medium, species)
    k_eff_lows = ωs./eff_medium.c

    k_eff_φs = wavenumber_low_volfrac(ωs, medium, species)
    k_effs = wavenumber(ωs, medium, species)

    @test norm(k_effs - k_eff_lows)/norm(k_effs) < 1e-5
    @test norm(k_effs[1] - k_eff_lows[1])/norm(k_effs[1]) < 1e-8
    @test norm(k_effs - k_eff_φs)/norm(k_effs) < 0.01
end

end
