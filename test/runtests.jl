using EffectiveWaves, Test
using LinearAlgebra, Statistics
using ClassicalOrthogonalPolynomials

@testset "Check material definitions" begin
    include("../src/acoustics/acoustic_mediums.jl")
    @test true
end

include("complex.jl")

# Single effective wavenumber tests
    include("acoustics-2D/strong_low_freq_effective.jl")
    include("acoustics-2D/high_frequency_effective.jl")

    include("acoustics-2D/large_vol_low_freq_effective.jl")
    include("acoustics-2D/weak_scatterers_effective.jl")

    include("acoustics-2D/pair-correlation.jl")
    include("acoustics-2D/average_integrand_kernel.jl")


# Eigensystems for different symmetries (i.e. a sphere fille with particles or a halfspace filled with particles)
    include("acoustics-3D/equivalent-symmetries.jl")
    include("acoustics-3D/low_volume_fraction.jl")
    include("acoustics-3D/low_frequency.jl")
    include("acoustics-3D/planar-symmetry.jl")
    include("acoustics-3D/sphere.jl")

    include("pair-correlation.jl")

# Eigensystems and fields for the case of two media
    include("acoustics-3D/Two-media/planar-symmetry.jl")

# test equivalence between methods for finding wavenumbers
    # test does not run on Julia version < 0.7 due to differences in Optim versions
    # include("path_mesh_wavenumbers.jl")

# Test functions used for the discretisation
    include("discretisation/numerical_integration.jl")

    # Test discretisation for 3D
    include("discretisation/function-approximations.jl")

# Test matching method
    include("match-wave/match_low_volumefraction.jl")
    include("match-wave/match_low_frequency.jl")
    include("match-wave/match_discrete_wave.jl")

# Test wiener hopf method for monopole scatterers
    # include("wiener-hopf-reflection.jl")
    # include("wiener-hopf-wave.jl")
