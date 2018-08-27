# Average waves are usually numerical solutions of the ensemble average wave equations. Here we mostly deal with the scenario where a halfspace filled with particles are ensemble averaged. The resulting govering equations are a 1D Fredholm integral equation.

export reflection_coefficient_integrated
export intergrand_kernel, average_wave_system, integrate_S, integrate_B
export integration_scheme, trap_scheme, simpson_scheme

include("average_waves.jl")
include("numerical_integration.jl")
include("integral_form.jl")
include("reflection.jl")
