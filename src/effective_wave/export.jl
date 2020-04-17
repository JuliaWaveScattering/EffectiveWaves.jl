# Effective waves are related to the ansatz u ~ amps.*exp(im*k_eff*x), where k_eff is the effective wavenumber with Im(k_eff) > 0.

# methods to find the wavenumbers.
export wavenumbers_bisection, wavenumbers_bisection_robust, wavenumbers_path, wavenumbers_mesh, reduce_kvecs, box_keff, NelderMeadparameters

export gray_square!, gray_square

export reflection_coefficient_low_volumefraction, wavenumber_low_volumefraction, wavenumber_very_low_volumefraction

include("wavemode.jl")
include("reflection_coefficient.jl")
include("dispersion.jl")
include("wavenumbers.jl")
include("solvers/wavenumbers_bisection.jl")
include("solvers/wavenumber_path.jl")
include("solvers/wavenumber_mesh.jl")
include("solvers/utils.jl")
