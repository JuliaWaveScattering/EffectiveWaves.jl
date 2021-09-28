# Average waves are usually numerical solutions of the ensemble average wave equations. Here we mostly deal with the scenario where a halfspace filled with particles are ensemble averaged. The resulting govering equations are a 1D Fredholm integral equation.

export intergrand_kernel, discrete_wave_system, discrete_system, integrate_S, integrate_B
export integration_scheme, trapezoidal_scheme, simpson_scheme
export discretewave_error, x_mesh

include("discrete_waves.jl")

include("planewaves/wave_types.jl")
include("planewaves/discrete_systems.jl")
include("planewaves/integral_form.jl")
include("planewaves/reflection.jl")

include("regular/discrete_systems.jl")

include("numerical_integration.jl")
