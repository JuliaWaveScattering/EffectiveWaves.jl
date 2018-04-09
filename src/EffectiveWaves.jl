# Loads all files
module EffectiveWaves

using BlackBoxOptim, Memoize

export  Specie, Medium, volume_fraction, Zn, p_speed,
        multispecies_wavenumber, wavenumber_very_low_volfrac, multispecies_challis, one_species_low_wavenumber,
        two_species_approx_wavenumber,
        opt_methods, optimal_species,
        gray_square!, gray_square,
        Brick, IronArmco, LeadAnnealed, RubberGum, FusedSilica, GlassPyrex,
        ClayRock, WaterDistilled, Glycerol, Hexadecane, Acetone, Benzene,
        Nitrobenzene, OliveOil, SodiumNitrate, AirDry,
        LimeStone, Clay, Calcite, SilicaQuartz


import Base.isequal, Base.(==)
import SpecialFunctions: besselj, hankelh1
try import BlackBoxOptim end

# using RecipesBase # Have not really needed yet

# push!(LOAD_PATH,"$(homedir())/.julia/v0.6/EffectiveWaves/examples/")

include("plot/graphics.jl")
include("particle.jl")
include("optimise_wavenumber.jl")
include("low_volfrac.jl")
include("wavenumber_challis.jl")
include("two_species_approximate.jl")
include("../examples/materials.jl")

end # module
