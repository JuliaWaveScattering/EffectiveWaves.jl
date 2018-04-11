import Base.Test: @testset, @test, @test_throws

using EffectiveWaves

# background medium
medium = Medium(1.0,1.0+0.0im)

# angular frequencies
ωs = [0.01,5.,20.,40.]
ωs2 = [60.,80.,120.]

@testset "Summary" begin

    @testset "Wavenumbers" begin

        @testset "weak scatterers" begin

        # Weak scatterers
        species = [
            Specie(ρ=10.,r=0.1, c=12., volfrac=0.1),
            Specie(ρ=3., r=0.2, c=2.0, volfrac=0.1)
        ]

        (β_eff,ρ_eff) = effective_material_properties(medium, species)
        k_eff_lows = ωs.*sqrt(ρ_eff/β_eff)
        k_eff_φs = wavenumber_low_volfrac(ωs, medium, species)
        k_effs = [wavenumber(ω, medium, species) for ω in ωs]

        @test norm(k_effs - k_eff_φs) < 0.002*norm(k_effs)
        @test norm(k_effs - k_eff_lows) < 0.06*norm(k_effs)
        @test norm(k_effs[1] - k_eff_lows[1]) < 0.02*norm(k_effs[1])

        # Large weak scatterers with low volume fraciton
        species = [
            Specie(ρ=10.,r=1.9, c=12., volfrac=0.04),
            Specie(ρ=3., r=0.7, c=2.0, volfrac=0.02)
        ]
        end

        @testset "high frequency" begin

        k_eff_φs = wavenumber_low_volfrac(ωs2, medium, species)
        k_effs = [wavenumber(ω, medium, species) for ω in ωs2]

        @test norm(k_effs - k_eff_φs) < 1e-4*norm(k_effs)

        end

        @testset "large volume fraction and low frequency" begin

        # large volume fraction scatterers,  small size amd on strong scatterer. This is a problamatic case.
        species = [
            Specie(ρ=5.,r=0.004, c=1.2, volfrac=0.4),
            Specie(ρ=0.3, r=0.002, c=0.4, volfrac=0.3)
        ]

        (β_eff,ρ_eff) = effective_material_properties(medium, species)
        k_eff_lows = ωs.*sqrt(ρ_eff/β_eff)
        k_effs = [wavenumber(ω, medium, species) for ω in ωs]

        @test norm(k_effs - k_eff_lows) < 0.02*norm(k_effs)
        @test norm(k_effs[1] - k_eff_lows[1]) < 0.01*norm(k_effs)

        end
    end

end
