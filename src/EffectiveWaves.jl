# Loads all files
module EffectiveWaves


export  Specie, Medium, volume_fraction, Zn, p_speed,
        multispecies_wavenumber, multispecies_wavenumber_low_volfrac, multispecies_challis, one_species_low_wavenumber,
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

include("plot/graphics.jl")
include("scattering.jl")
include("multi-species.jl")
include("multi-species_challis.jl")
include("two_species_approximate.jl")
include("../examples/materials.jl")

end # module
