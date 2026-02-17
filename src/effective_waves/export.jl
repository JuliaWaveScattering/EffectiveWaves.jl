# Effective waves are related to the ansatz u ~ amps.*exp(im*k_eff*x), where k_eff is the effective wavenumber with Im(k_eff) > 0.

# methods to find the wavenumbers.
export wavenumbers_bisection, wavenumbers_bisection_robust, wavenumbers_path, wavenumbers_mesh, reduce_kvecs, box_keff, NelderMeadparameters

export gray_square!, gray_square

export reflection_coefficient_low_volumefraction, wavenumber_low_volumefraction, wavenumber_very_low_volumefraction

# include("effective_wave/planewaves/eigensystem.jl")
include("plane_waves/plane-eigensystem.jl")
# include("effective_wave/regular/eigensystem.jl")
include("plane_waves/boundary-condition.jl")
include("plane_waves/reflection-transmission.jl")
include("plane_waves/reflection_coefficient.jl")

include("regular_waves/regular-eigensystem.jl")
include("regular_waves/boundary_condition.jl")
include("regular_waves/material_scattering_coefficients.jl")


include("dispersion.jl")
include("wavemode.jl")
include("asymptotics.jl")
include("wavenumbers.jl")

include("solvers/wavenumbers_bisection.jl")
include("solvers/wavenumbers_bisection_robust.jl")
include("solvers/wavenumber_path.jl")
include("solvers/utils.jl")
