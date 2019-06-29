# Loads all files
module EffectiveWaves

export  EffectiveWave, AverageWave # the two main types
export  MatchWave # a combination of the two types above

# for MatchWave
export  match_error, x_mesh_match

# for discrete method
export  x_mesh

# spherical bessel and hankel functions
export  sbesselj, shankelh1, diffsbessel, diffbessel

# for material and particle properties
export  Specie, Medium, volume_fraction, Zn, t_vectors, Nn, p_speed, maximum_hankel_order

# for effective waves
export  wavenumbers, wavenumber, effective_waves, transmission_angle,
        effective_wavevectors, scattering_amplitudes_average, scale_amplitudes_effective

export  reflection_coefficient, reflection_coefficients

export  effective_medium, reflection_coefficient_halfspace

# List of shorthand for some materials
export  Brick, IronArmco, LeadAnnealed, RubberGum, FusedSilica, GlassPyrex,
        ClayRock, WaterDistilled, Glycerol, Hexadecane, Acetone, Benzene,
        Nitrobenzene, OliveOil, SodiumNitrate, AirDry,
        LimeStone, Clay, Calcite, SilicaQuartz

import Base.isequal, Base.(==), Base.zero
import SpecialFunctions: besselj, hankelh1
import Statistics: mean, std

# try using BlackBoxOptim end
# using IterTools
using RecipesBase
using Optim
using OffsetArrays
using ApproxFun
using LinearAlgebra

# using RecipesBase # Have not really needed yet

# push!(LOAD_PATH,"$(homedir())/.julia/v0.6/EffectiveWaves/examples/")

include("specialfunctions.jl")
include("particle.jl")
include("t-matrix.jl")

include("effective_waves/effective_waves_export.jl")
include("average_waves/average_waves_export.jl")
include("match_waves/match_waves.jl")
include("match_waves/match_arrays.jl")
include("match_waves/reflection.jl")

include("../examples/materials.jl")

include("plot/graphics.jl")
include("plot/plot.jl")

end # module
