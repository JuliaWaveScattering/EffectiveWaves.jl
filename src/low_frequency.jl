"the effective low frequency bulk modulus and density of a material filled with particles"
function effective_material_properties{T<:Number}(medium::Medium{T}, species::Array{Specie{T}})
    # total number density
    η = sum(s.num_density for s in species)

    φ = sum(volume_fraction.(species))

    # calculate effective properties
    β = medium.ρ*medium.c^2 # medium bulk modulus
    β_eff = real( (1-φ)/β + (φ/η)*sum(s.num_density/(s.ρ*s.c^2) for s in species) )^(-1.0)
    ρ_frac = (φ/η)*sum(s.num_density*(medium.ρ - s.ρ)/(medium.ρ + s.ρ) for s in species)
    ρ_eff = medium.ρ*(1 - ρ_frac)/(1 + ρ_frac)

    return (β_eff, ρ_eff)
end

"the reflection from a homogenious halfspace, which is also the low frequency reflection from a particulate material when using the effective_material_properties."
function reflection_coefficient_halfspace{T<:Number}(incident_medium::Medium{T}, reflect_medium::Medium{T}; θ_inc::T = zero(T), tol = 1e-6)
    q = real(reflect_medium.c*reflect_medium.ρ/(incident_medium.c*incident_medium.ρ))
    snell(θ) = abs(reflect_medium.c*sin(θ_inc) - incident_medium.c*sin(θ))
    result = optimize(snell, [θ_inc,0.]; g_tol = tol^2.0)
    θ_trans = result.minimizer[1] + im*result.minimizer[2]
    R = (q*cos(θ_inc) - cos(θ_trans))/(q*cos(θ_inc) + cos(θ_trans))
end
