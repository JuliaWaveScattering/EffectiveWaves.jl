using EffectiveWaves, Test
using LinearAlgebra

@testset "Check material definitions" begin
    include("../src/materials.jl")
    @test true
end

@testset "Examples from: Reflection from multi-species.., Proc.R.Soc.(2018)" begin

    include("../docs/src/examples/concrete/concrete_species.jl")
    include("../docs/src/examples/concrete/concrete_species_volfrac.jl")
    # include("../docs/src/examples/concrete/concrete_species_large-freq.jl") # takes longer

    # Takes too long
    include("../docs/src/examples/emulsion/fluid_species.jl")
    include("../docs/src/examples/emulsion/fluid_species_volfrac.jl")
    # include("../docs/src/examples/emulsion/fluid_species_large-freq.jl") # takes longer

    @test true
end

include("specialfunctions.jl")

include("types_constructors.jl")

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
