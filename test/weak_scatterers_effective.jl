@testset "weak scatterers" begin
    # background medium
    medium = Medium(1.0,1.0+0.0im)

    # angular frequencies
    ωs = [0.001,20.]

    species = [
        Specie(ρ=10.,r=0.01, c=12., volfrac=0.05),
        Specie(ρ=3., r=0.2, c=2.0, volfrac=0.04)
    ]

    # wavenumbers
    eff_medium = effective_medium(medium, species)

    k_eff_lows = ωs./eff_medium.c
    k_eff_φs = wavenumber_low_volfrac(ωs, medium, species; tol=1e-9)

    k_effs = [
        wavenumbers(ω, medium, species;
            num_wavenumbers=1, tol = 1e-8)
    for ω in ωs]
    inds = [argmin(abs.(k_effs[i] .- k_eff_φs[i])) for i in eachindex(ωs)]
    k_effs2 = [k_effs[i][inds[i]] for i in eachindex(inds)]

    @test norm(k_effs2 - k_eff_φs)/norm(k_eff_φs) < 0.002
    @test norm(k_effs2 - k_eff_lows)/norm(k_eff_lows) < 0.01
    @test norm(k_effs2[1] - k_eff_lows[1])/norm(k_effs[1]) < 4e-7

    wave_effs_2 = [EffectiveWave(ωs[i], k_effs2[i], medium, species; tol = 1e-9) for i in eachindex(ωs)]
    wave_effs_φs = [EffectiveWave(ωs[i], k_eff_φs[i], medium, species; tol = 1e-9) for i in eachindex(ωs)]

    # reflection coefficient
    Rs2 = reflection_coefficients(ωs, medium, species; tol=1e-9, num_wavenumbers=1)

    # if the wavenumber k_effs are already calculated, the below is faster
    # pairs each ω in ωs with each wave in wave_effs2 to calculate the refleciton coefficients

    Rs = map(eachindex(ωs)) do i
        w_effs = [EffectiveWave(ωs[i], k_eff, medium, species; tol = 1e-9) for k_eff in k_effs[i]]
        reflection_coefficient(ωs[i], w_effs, medium, species; tol = 1e-9)
    end

    @test norm(Rs - Rs2)/norm(Rs) < 5e-6

    Rs1 = map(eachindex(ωs)) do i
        w_effs = EffectiveWave(ωs[i], k_effs[i][1], medium, species; tol = 1e-9)
        reflection_coefficient(ωs[i], w_effs, medium, species)
    end

    # Direct incidence
    R_low = reflection_coefficient_halfspace(medium, eff_medium; tol=1e-9)
    Rs_φs = reflection_coefficients(ωs, wave_effs_φs, medium, species; tol=1e-9)
    # the below takes a low-volfrac expansion for both the wavenumber and reflection coefficient
    Rs_φs2 = reflection_coefficient_low_volfrac(ωs, medium, species; tol=1e-9)

    len = length(ωs)
    @test maximum(abs.(Rs_φs - Rs1)) < 2e-4 # already relative to incident wave amplitude = 1
    @test maximum(abs.(Rs_φs2 - Rs1)) < 1e-3
    @test abs(R_low - Rs1[1]) < 1e-7

    # Vary angle of incidence θin
    θs = 0.1:0.3:(π/2)
    # θs = (-π/1.9):0.3:(π/1.9)
    R_low = [reflection_coefficient_halfspace(medium, eff_medium; θin = θ) for θ in θs]

    Rs = map(θs) do θ  # , hankel_order =7
        wave_effs_2 = [EffectiveWave(ωs[i], k_effs2[i], medium, species; tol = 1e-9, θin = θ) for i in eachindex(ωs)]
        reflection_coefficients(ωs, wave_effs_2, medium, species; θin = θ)
    end
    Rs_φs = [reflection_coefficient_low_volfrac(ωs, medium, species; θin = θ, tol = 1e-9, hankel_order =7) for θ in θs];

    @test maximum(abs(R_low[i] - Rs[i][1]) for i in 1:length(R_low)) < 2e-6
    @test maximum(norm(R)/len for R in (Rs_φs - Rs)) < 0.01
end
