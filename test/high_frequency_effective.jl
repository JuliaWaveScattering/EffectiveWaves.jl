using EffectiveWaves, Test
using LinearAlgebra

@testset "high frequency effective" begin
        medium = Acoustic(2; ρ=1.0, c=1.0)

        # Large weak scatterers with low volume fraciton
        ms = MultipleScattering

        p1 = Particle(Acoustic(2; ρ=10.0, c=12.0),ms.Circle(1.9))
        p2 = Particle(Acoustic(2; ρ=3.0, c=2.0),ms.Circle(0.7))

        species = [
            Specie(p1; volume_fraction=0.04),
            Specie(p2; volume_fraction=0.02)
        ]
        ωs2 = [120.]

        tol = 1e-6
        k_eff_φs = wavenumber_low_volumefraction(ωs2, medium, species; tol=tol)
        k_effs = [wavenumbers(ω, medium, species; tol=tol, num_wavenumbers=1) for ω in ωs2]
        inds = [argmin(abs.(k_effs[i] .- k_eff_φs[i])) for i in eachindex(ωs2)]
        k_effs2 = [k_effs[i][inds[i]] for i in eachindex(inds)]

        @test norm(k_effs2 - k_eff_φs)/norm(k_effs2) < tol

        Rs = map(eachindex(ωs2)) do i
            wave = EffectivePlaneWaveMode(ωs2[i], k_effs2[i], medium, species)
            reflection_coefficient(ωs2[i], wave, medium, species)
        end
        # warning is expected, as k_eff_φs are assymptotic approximations.
        Rs_φs = map(eachindex(ωs2)) do i
            wave = EffectivePlaneWaveMode(ωs2[i], k_eff_φs[i], medium, species)
            reflection_coefficient(ωs2[i], wave, medium, species)
        end
        Rs_φs2 = reflection_coefficient_low_volumefraction(ωs2, medium, species)

        # the incident wave has amplitude 1, so this is already a relative difference
        @test maximum(abs.(Rs_φs - Rs)) < tol
        @test maximum(abs.(Rs_φs2 - Rs)) < tol
end
