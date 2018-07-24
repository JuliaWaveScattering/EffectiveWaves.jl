@testset "high frequency effective" begin
    medium = Medium(1.0,1.0+0.0im)
    # Large weak scatterers with low volume fraciton
    species = [
        Specie(ρ=10.,r=1.9, c=12., volfrac=0.04),
        Specie(ρ=3., r=0.7, c=2.0, volfrac=0.02)
    ]
    ωs2 = 20.:30.:121

    k_eff_φs = wavenumber_low_volfrac(ωs2, medium, species; tol=1e-5) # lowering tol to speed up calculation
    k_effs = [wavenumbers(ω, medium, species; tol=1e-6) for ω in ωs2]
    # k_effs = wavenumber(ωs2, medium, species; tol=1e-7) # lowering tol to speed up calculation
    inds = [indmin(abs.(k)) for k in (k_effs .- k_eff_φs)]
    k_effs2 = [k_effs[i][inds[i]] for i in eachindex(inds)]

    @test norm(k_effs2 - k_eff_φs)/norm(k_effs) < 2e-4

    Rs = map(eachindex(ωs2)) do i
        wave = EffectiveWave(ωs2[i], k_effs2[i], medium, species)
        reflection_coefficient(ωs2[i], wave, medium, species)
    end
    Rs_φs = map(eachindex(ωs2)) do i
        wave = EffectiveWave(ωs2[i], k_eff_φs[i], medium, species)
        reflection_coefficient(ωs2[i], wave, medium, species)
    end
    Rs_φs2 = reflection_coefficient_low_volfrac(ωs2, medium, species)

    len = length(ωs2)
    @test norm(Rs_φs - Rs)/len < 5e-5
    @test norm(Rs_φs2 - Rs)/len < 5e-5 # the incident wave has amplitude 1, so this is a tiny difference
end
