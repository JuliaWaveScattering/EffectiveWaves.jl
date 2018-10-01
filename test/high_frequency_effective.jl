@testset "high frequency effective" begin
    medium = Medium(1.0,1.0+0.0im)
    # Large weak scatterers with low volume fraciton
    species = [
        Specie(ρ=10.,r=1.9, c=12., volfrac=0.04),
        Specie(ρ=3., r=0.7, c=2.0, volfrac=0.02)
    ]
    ωs2 = [120.]

    tol = 1e-6
    k_eff_φs = wavenumber_low_volfrac(ωs2, medium, species; tol=tol)
    k_effs = [wavenumbers(ω, medium, species; tol=tol, num_wavenumbers=1) for ω in ωs2]
    inds = [indmin(abs.(k)) for k in (k_effs .- k_eff_φs)]
    k_effs2 = [k_effs[i][inds[i]] for i in eachindex(inds)]

    @test norm(k_effs2 - k_eff_φs)/norm(k_effs) < tol

    Rs = map(eachindex(ωs2)) do i
        wave = EffectiveWave(ωs2[i], k_effs2[i], medium, species)
        reflection_coefficient(ωs2[i], wave, medium, species)
    end
    # warning is expected, as k_eff_φs are assymptotic approximations.
    Rs_φs = map(eachindex(ωs2)) do i
        wave = EffectiveWave(ωs2[i], k_eff_φs[i], medium, species)
        reflection_coefficient(ωs2[i], wave, medium, species)
    end
    Rs_φs2 = reflection_coefficient_low_volfrac(ωs2, medium, species)

    @test maximum(abs.(Rs_φs - Rs)) < tol
    @test maximum(abs.(Rs_φs - Rs)) < tol # the incident wave has amplitude 1, so this is already a relative difference
end
