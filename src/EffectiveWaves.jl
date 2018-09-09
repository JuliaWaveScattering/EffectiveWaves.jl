# Loads all files
module EffectiveWaves

export  EffectiveWave, AverageWave # the two main types
export  MatchWave # a combination of the two types above

export  Specie, Medium, volume_fraction, Zn, t_vectors, Nn, p_speed, maximum_hankel_order

export  wavenumbers, wavenumber, effective_waves, transmission_angle,
        reduced_amplitudes_effective, scattering_amplitudes_average, scale_amplitudes_effective

export  reflection_coefficient, reflection_coefficients

export  effective_medium, reflection_coefficient_halfspace

export  x_mesh, x_mesh_match

export  Brick, IronArmco, LeadAnnealed, RubberGum, FusedSilica, GlassPyrex,
        ClayRock, WaterDistilled, Glycerol, Hexadecane, Acetone, Benzene,
        Nitrobenzene, OliveOil, SodiumNitrate, AirDry,
        LimeStone, Clay, Calcite, SilicaQuartz

import Base.isequal, Base.(==), Base.zero
import SpecialFunctions: besselj, hankelh1

# try using BlackBoxOptim end
using RecipesBase
using Optim
using NLsolve
# using IterTools
using OffsetArrays
using ApproxFun

# using RecipesBase # Have not really needed yet

# push!(LOAD_PATH,"$(homedir())/.julia/v0.6/EffectiveWaves/examples/")

include("particle.jl")

include("effective_waves/effective_waves_export.jl")
include("average_waves/average_waves_export.jl")
include("match_waves/match_waves.jl")
include("match_waves/match_arrays.jl")

include("../examples/materials.jl")

include("plot/graphics.jl")
include("plot/plot.jl")

end # module
