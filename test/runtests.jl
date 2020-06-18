using EffectiveWaves, Test
using LinearAlgebra

@testset "Check material definitions" begin
    include("../src/acoustics/acoustic_mediums.jl")
    @test true
end

include("complex.jl")

# Single effective wavenumber tests
    @time include("effective-2D/strong_low_freq_effective.jl")
    @time include("effective-2D/high_frequency_effective.jl")

    @time include("effective-2D/large_vol_low_freq_effective.jl")
    @time include("effective-2D/weak_scatterers_effective.jl")


# Eigensystems for different symmetries (i.e. a sphere fille with particles or a halfspace filled with particles)
    @time include("effective-3D/equivalent-symmetries.jl")
    @time include("effective-3D/fourth-order-matrix.jl")
    @time include("effective-3D/low_volume_fraction.jl")
    @time include("effective-3D/low_frequency.jl")

# test equivalence between methods for finding wavenumbers
    # test does not run on Julia version < 0.7 due to differences in Optim versions
    include("path_mesh_wavenumbers.jl")

# Test functions used for the discretisation
    include("discretisation/numerical_integration.jl")
    include("discretisation/integrated_reflection.jl")
    include("discretisation/average_integrand_kernel.jl")

# Test matching method
    include("match-wave/match_low_volumefraction.jl")
    include("match-wave/match_low_frequency.jl")
    include("match-wave/match_discrete_wave.jl")

# Test wiener hopf method for monopole scatterers
    include("wiener-hopf-reflection.jl")
    include("wiener-hopf-wave.jl")
