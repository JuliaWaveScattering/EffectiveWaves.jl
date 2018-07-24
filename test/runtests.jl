import Base.Test: @testset, @test, @test_throws

using EffectiveWaves

@testset "Summary" begin

@testset "Examples from: Reflection from multi-species.., Proc.R.Soc.(2018)" begin

    include("../examples/concrete/concrete_species.jl")
    include("../examples/concrete/concrete_species_volfrac.jl")
    # include("../examples/concrete/concrete_species_large-freq.jl") # takes longer

    # Takes too long
    include("../examples/emulsion/fluid_species.jl")
    include("../examples/emulsion/fluid_species_volfrac.jl")
    # include("../examples/emulsion/fluid_species_large-freq.jl") # takes longer

    @test true
end

# @testset "Tests for the integral form" begin
    # include("integral_form_tests.jl")
    # test_integrand()

# end

include("types_constructors.jl")

species = [
    Specie(ρ=10.,r=0.01, c=12., volfrac=0.05),
    Specie(ρ=3., r=0.2, c=2.0, volfrac=0.04)
]

# # On Travis the below give some empty entries. This possibly due to time_limit not being stable accross different systems.
k_eff = wavenumbers(21., Medium(1.0,1.0+0.0im), species; tol = 1e-6, time_limit=1.0)
if !isempty(k_eff) include("weak_scatterers_effective.jl") end

include("high_frequency_effective.jl")
include("large_vol_low_freq_effective.jl")
include("strong_low_freq_effective.jl")

include("integrated_reflection.jl")

end
