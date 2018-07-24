# Loads all files
module EffectiveWaves

export  EffectiveWave, AverageWave # the two main types

export  Specie, Medium, volume_fraction, Zn, Nn, p_speed, maximum_hankel_order

export  wavenumbers, wavenumber, reflection_coefficient, transmission_angle,
        reduced_amplitudes_effective, scattering_amplitudes_average, scale_amplitudes_effective

export  effective_medium, reflection_coefficient_halfspace

export  Brick, IronArmco, LeadAnnealed, RubberGum, FusedSilica, GlassPyrex,
        ClayRock, WaterDistilled, Glycerol, Hexadecane, Acetone, Benzene,
        Nitrobenzene, OliveOil, SodiumNitrate, AirDry,
        LimeStone, Clay, Calcite, SilicaQuartz

import Base.isequal, Base.(==)
import SpecialFunctions: besselj, hankelh1

# try using BlackBoxOptim end

using Optim
using IterTools
using OffsetArrays
using ApproxFun

# using RecipesBase # Have not really needed yet

# push!(LOAD_PATH,"$(homedir())/.julia/v0.6/EffectiveWaves/examples/")

include("plot/graphics.jl")
include("particle.jl")

include("effective_waves/effective_waves_export.jl")

include("average_waves/average_waves_export.jl")

include("../examples/materials.jl")

end # module
