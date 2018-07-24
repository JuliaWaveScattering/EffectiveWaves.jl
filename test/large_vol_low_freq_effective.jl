# large volume fraction scatterers,  small size amd on strong scatterer. This is a problamatic case.
@testset "large volume fraction and low frequency effective" begin
    medium = Medium(1.0,1.0+0.0im)
    species = [
        Specie(ρ=5.,r=0.004, c=1.2, volfrac=0.4),
        Specie(ρ=0.3, r=0.002, c=0.4, volfrac=0.3)
    ]
    # angular frequencies
    ωs = 0.001:5.0:21.

    eff_medium = effective_medium(medium, species)
    k_eff_lows = ωs./eff_medium.c
    # k_effs = wavenumber(ωs, medium, species)
    k_effs = [wavenumbers(ω, medium, species) for ω in ωs]
    inds = [indmin(abs.(k)) for k in (k_effs .- k_eff_lows)]
    k_effs2 = [k_effs[i][inds[i]] for i in eachindex(inds)]

    @test norm(k_effs2 .- k_eff_lows)/norm(k_eff_lows) < 0.01
    @test norm(k_effs2[1] - k_eff_lows[1])/norm(k_eff_lows[1]) < 2e-6

    Rs = map(eachindex(ωs)) do i
        wave = EffectiveWave(ωs[i], k_effs2[i], medium, species)
        reflection_coefficient(ωs[i], wave, medium, species)
    end
    R_low2 = map(eachindex(ωs)) do i
        wave = EffectiveWave(ωs[i], k_eff_lows[i], medium, species)
        reflection_coefficient(ωs[i], wave, medium, species)
    end
    R_low = reflection_coefficient_halfspace(medium, eff_medium)

    @test norm(R_low - Rs[1]) < 5e-7
    @test norm(R_low2[1] - Rs[1]) < 5e-7
    @test norm(R_low2 - Rs) < 0.005
end
