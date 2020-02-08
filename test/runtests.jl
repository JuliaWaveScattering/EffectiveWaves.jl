using EffectiveWaves, Test
using LinearAlgebra

@testset "Check material definitions" begin
    include("../src/acoustics/acoustic_mediums.jl")
    @test true
end

include("specialfunctions.jl")
include("complex.jl")

# Single effective wavenumber tests
@time    include("strong_low_freq_effective.jl")
    @time include("high_frequency_effective.jl")

    @time include("large_vol_low_freq_effective.jl")
    @time include("weak_scatterers_effective.jl")

# test equivalence between methods for finding wavenumbers
    # test does not run on Julia version < 0.7 due to differences in Optim versions
    include("path_mesh_wavenumbers.jl")

# Test functions used for the discretisation
    include("numerical_integration.jl")
    include("integrated_reflection.jl")
    include("average_integrand_kernel.jl")

# Test matching method
    include("match_wave.jl")

# Test wiener hopf method for monopole scatterers
    include("wiener-hopf-reflection.jl")
    include("wiener-hopf-wave.jl")
