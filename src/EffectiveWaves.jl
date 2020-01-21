# Loads all files
module EffectiveWaves

include("export_import.jl")
# try using BlackBoxOptim end
# using IterTools
using RecipesBase, OffsetArrays, LinearAlgebra
using Optim, ApproxFun # Heavy packages
using WignerSymbols, GSL

using Reexport
@reexport using MultipleScattering

include("specialfunctions.jl")

include("effective_wave/export.jl")
include("discrete_wave/export.jl")
include("match_waves/match_waves.jl")
include("match_waves/match_arrays.jl")
include("match_waves/reflection.jl")

include("acoustics/export.jl")

include("materials.jl")

include("plot/graphics.jl")
include("plot/plot.jl")

end # module
