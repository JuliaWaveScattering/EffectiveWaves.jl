# Effective waves are related to the ansatz u ~ amps.*exp(im*k_eff*x), where k_eff is the effective wavenumber with Im(k_eff) > 0.

export effectivewave_system # supplies a matrix used for the disperision equation and effective eignvectors
export dispersion_equation
export wavenumbers_path, wavenumbers_mesh, reduce_kvecs

export gray_square!, gray_square

export reflection_coefficient_low_volumefraction, wavenumber_low_volumefraction, wavenumber_very_low_volumefraction

include("effective_wave.jl")
include("wavemode.jl")
include("reflection_coefficient.jl")
include("dispersion.jl")
include("wavenumber_effective.jl")
include("wavenumber_path.jl")
include("wavenumber_mesh.jl")
include("utils.jl")
