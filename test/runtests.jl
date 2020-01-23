using EffectiveWaves, Test
using LinearAlgebra

@testset "Check material definitions" begin
    include("../src/acoustics/acoustic_mediums.jl")
    @test true
end

include("specialfunctions.jl")
include("complex.jl")

@testset "Single effective wave" begin
    include("strong_low_freq_effective.jl")
    include("high_frequency_effective.jl")

    include("large_vol_low_freq_effective.jl")
    include("weak_scatterers_effective.jl")
end

# test does not run on Julia version < 0.7 due to differences in Optim versions
include("path_mesh_wavenumbers.jl")

include("numerical_integration.jl")
include("integrated_reflection.jl")
include("average_integrand_kernel.jl")

include("match_wave.jl")
include("wiener-hopf-reflection.jl")
include("wiener-hopf-wave.jl")
