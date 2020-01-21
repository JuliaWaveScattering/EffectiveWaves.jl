# Here are the main exported functions and types. Note there are other exported functions and types in files such as "effective_waves/effective_waves_export" and "discrete_wave/export.jl."

export  EffectivePlaneWaveMode, DiscretePlaneWaveMode # the two main types
export  MatchPlaneWaveMode # a combination of the two types above

export Specie, Species, number_density, volume_fraction
export Material

# for MatchPlaneWaveMode
export  match_error, x_mesh_match

# for discrete method
export  x_mesh

# for material and particle properties
# export  Specie, Medium, volume_fraction, Zn, t_vectors, Nn, p_speed, maximum_basis_order
export  Nn, p_speed#, maximum_basis_order

# for effective waves
export  wavenumbers, wavenumber, effective_waves, transmission_angle,
        effective_wavemodes, scattering_amplitudes_average, scale_amplitudes_effective

export  reflection_coefficient, reflection_coefficients

export  effective_medium

# List of shorthand for some materials
export  Brick, IronArmco, LeadAnnealed, RubberGum, FusedSilica, GlassPyrex,
        ClayRock, WaterDistilled, Glycerol, Hexadecane, Acetone, Benzene,
        Nitrobenzene, OliveOil, SodiumNitrate, AirDry,
        LimeStone, Clay, Calcite, SilicaQuartz

import Base.isequal, Base.(==), Base.zero
import SpecialFunctions: besselj, hankelh1
import Statistics: mean, std
