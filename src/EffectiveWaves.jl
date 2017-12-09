# Loads all files
module EffectiveWaves


export  Specie, Medium, volume_fraction, Zn, p_speed
        multispecies_wavenumber, multispecies_wavenumber_low_volfrac, multispecies_challis, one_species_low_wavenumber,
        opt_methods, optimal_species,
        gray_square!, gray_square


import Base.isequal, Base.(==)
import SpecialFunctions: besselj, hankelh1
try import BlackBoxOptim end

include("scattering.jl")
include("multi-species.jl")
include("multi-species_challis.jl")
include("two_species_approximate.jl")
include("plot/graphics.jl")
include("../examples/materials.jl")

end # module
