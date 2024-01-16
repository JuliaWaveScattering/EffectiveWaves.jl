export p_speed
export far_field_pattern, pair_field_pattern, diff_far_field_pattern
export wavenumber_challis, one_species_low_wavenumber,
        two_species_approx_wavenumber
export wienerhopf_reflection_coefficient, wienerhopf_mode_amplitudes

include("acoustic_mediums.jl")

include("low_frequency.jl")

include("effective_wave/special_functions.jl")
# include("effective_wave/asymptotics.jl")

# include("effective_wave/planewaves/boundary_condition.jl")
# include("effective_wave/planewaves/wavemode-wienerhopf.jl")

# regular means a smooth field assumption and includes cases such as a sphere filled with particles.

include("effective_wave/far_fields.jl")
include("effective_wave/low_volumefraction.jl")
include("effective_wave/alternative_wavenumbers.jl")
include("effective_wave/two_species_approximate.jl")


include("discrete_wave/export.jl")
