# Loads all files
module EffectiveWaves

include("export_import.jl")
# try using BlackBoxOptim end
# using IterTools
using RecipesBase
using Optim
using OffsetArrays
using ApproxFun
using LinearAlgebra
using WignerSymbols
using GSL

include("specialfunctions.jl")
include("particle.jl")
include("t-matrix.jl")

include("effective_waves/effective_waves_export.jl")
include("average_waves/average_waves_export.jl")
include("match_waves/match_waves.jl")
include("match_waves/match_arrays.jl")
include("match_waves/reflection.jl")

include("materials.jl")

include("plot/graphics.jl")
include("plot/plot.jl")

end # module
