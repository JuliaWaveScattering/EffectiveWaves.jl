
# Average waves are usually numerical solutions of the ensemble average wave equations. Here we mostly deal with the scenario where a halfspace filled with particles are ensemble averaged. The resulting govering equations are a 1D Fredholm integral equation.

export reflection_coefficient_integrated
export  trap_scheme, simpson_scheme #,intergrand_kernel, integral_form

include("average_waves.jl")
include("numerical_integration.jl")

include("reflection.jl")
