# Loads all files
module EffectiveWaves

# Here are the main exported functions and types. Note there are other exported functions and types in files such as "effective_wave/export.jl" and "discrete_wave/export.jl"

export  EffectivePlaneWaveMode, DiscretePlaneWaveMode # the two main types
export  MatchPlaneWaveMode # a combination of the two types above
export  Material, Specie, Species, SetupSymmetry, number_density, volume_fraction

# for MatchPlaneWaveMode
export  match_error, x_mesh_match
# export  Nn, p_speed#, maximum_basis_order

# for effective waves
export wavenumbers, wavenumber, eigenvectors
export wavemodes, wavemode, wavemode_wienerhopf
export solve_boundary_condition

export dispersion_equation, dispersion_complex, eigensystem # supplies a matrix used for the disperision equation and effective eignvectors

export scattering_amplitudes_average

export  reflection_coefficient, reflection_coefficients
export  effective_medium

# List of shorthand for some materials
export  Brick, IronArmco, LeadAnnealed, RubberGum, FusedSilica, GlassPyrex,
        ClayRock, WaterDistilled, Glycerol, Hexadecane, Acetone, Benzene,
        Nitrobenzene, OliveOil, SodiumNitrate, AirDry,
        LimeStone, Clay, Calcite, SilicaQuartz

import Base.isequal, Base.(==), Base.zero
import SpecialFunctions: besselj, hankelh1
import StaticArrays: SVector
import Statistics: mean, std

using Reexport
@reexport using MultipleScattering

using RecipesBase, OffsetArrays, LinearAlgebra

# Heavy package
using Optim: optimize, Optim, FixedParameters, Options, LBFGS, Fminbox, NelderMead
# Below doesn't work on Julia 1.1
# @reexport using Optim: optimize, FixedParameters, Options, LBFGS, Fminbox, NelderMead
import Optim: simplexer

# Heavy package
using ApproxFun: ApproxFun.(..), Fun, Segment, Domain, Chebyshev, DefiniteIntegral, LowRankFun, Interval, Legendre

using MDBM

include("material_types.jl")
include("wave_types.jl")
include("specialfunctions.jl")

include("effective_wave/export.jl")
include("acoustics/export.jl")

include("match_waves/match_waves.jl")
include("match_waves/match_arrays.jl")
include("match_waves/reflection.jl")

include("plot/graphics.jl")
include("plot/plot.jl")

end # module
