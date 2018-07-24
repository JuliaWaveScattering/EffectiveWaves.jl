import Base.Test: @testset, @test, @test_throws

using EffectiveWaves

@testset "Summary" begin

include("types_constructors.jl")

include("large_vol_low_freq_effective.jl")
include("integrated_reflection.jl")

end
