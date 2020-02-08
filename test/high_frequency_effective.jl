using EffectiveWaves, Test
using LinearAlgebra

@testset "high frequency effective" begin

        medium = Acoustic(2; ρ=1.0, c=1.0)
        basis_order = 3 # high frequency should use higher basis_order, but takes a while longer.

        # Large weak scatterers with low volume fraciton
        ms = MultipleScattering

        s1 = Specie(Acoustic(2; ρ=10.0, c=12.0),ms.Circle(1.9); volume_fraction=0.04)
        s2 = Specie(Acoustic(2; ρ=3.0, c=2.0),ms.Circle(0.7); volume_fraction=0.02)

        species = [s1, s2]

        ωs2 = [120.]

        tol = 1e-5 # low tolerance to speed up wavenumbers
        k_eff_φs = wavenumber_low_volumefraction(ωs2, medium, species; basis_order=basis_order)
        k_effs = [wavenumbers(ω, medium, species; tol=tol, num_wavenumbers=1, basis_order=basis_order) for ω in ωs2]

        inds = [argmin(abs.(k_effs[i] .- k_eff_φs[i])) for i in eachindex(ωs2)]
        k_effs2 = [k_effs[i][inds[i]] for i in eachindex(inds)]

        @test norm(k_effs2 - k_eff_φs)/norm(k_effs2) < tol

        # let's now assume the particles are fill a halfspace. Then we can calculate a reflection coefficent from this halfspace
        normal = [-1.0,0.0] # an outward normal to the surface
        material = Material(Halfspace(normal),species)

        # define a plane wave source travelling directly towards the material
        source = PlaneSource(medium, -normal)

        Rs = map(eachindex(ωs2)) do i
            wave = effective_wavemode(ωs2[i], k_effs2[i], source, material; basis_order=basis_order)
            reflection_coefficient(ωs2[i], wave, source, material)
        end
        # warning is expected, as k_eff_φs are assymptotic approximations.
        Rs_φs = map(eachindex(ωs2)) do i
            wave = effective_wavemode(ωs2[i], k_eff_φs[i], source, material; basis_order=basis_order)
            reflection_coefficient(ωs2[i], wave, source, material)
        end
        Rs_φs2 = reflection_coefficient_low_volumefraction(ωs2, source, material; basis_order=basis_order)

        # the incident wave has amplitude 1, so this is already a relative difference
        @test maximum(abs.(Rs_φs - Rs)) < 1e-10
        @test maximum(abs.(Rs_φs2 - Rs)) < 1e-8
end
