# Loads all files
module EffectiveWaves

export  Specie, Medium, volume_fraction, Zn, p_speed,
        wavenumber, reflection_coefficient, transmission_angle, transmission_scattering_coefficients,
        effective_medium, reflection_coefficient_halfspace,
        wavenumber_low_volfrac, wavenumber_very_low_volfrac, wavenumber_challis, one_species_low_wavenumber,
        two_species_approx_wavenumber,
        opt_methods, optimal_species,
        gray_square!, gray_square,
        Brick, IronArmco, LeadAnnealed, RubberGum, FusedSilica, GlassPyrex,
        ClayRock, WaterDistilled, Glycerol, Hexadecane, Acetone, Benzene,
        Nitrobenzene, OliveOil, SodiumNitrate, AirDry,
        LimeStone, Clay, Calcite, SilicaQuartz

import Base.isequal, Base.(==)
import SpecialFunctions: besselj, hankelh1

try using BlackBoxOptim end
using Memoize
using Optim

# using RecipesBase # Have not really needed yet

# push!(LOAD_PATH,"$(homedir())/.julia/v0.6/EffectiveWaves/examples/")

include("plot/graphics.jl")
include("particle.jl")
include("optimise_wavenumber.jl")
include("wavenumbers.jl")
include("reflection_transmission.jl")
include("low_frequency.jl")
include("low_volfrac.jl")
include("alternative_wavenumbers.jl")
include("two_species_approximate.jl")
include("../examples/materials.jl")

end # module
