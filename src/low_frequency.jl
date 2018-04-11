"the effective low frequency bulk modulus and density of a material filled with particles"
function effective_medium{T<:Number}(medium::Medium{T}, species::Array{Specie{T}})
    # total number density
    η = sum(s.num_density for s in species)

    φ = sum(volume_fraction.(species))

    # calculate effective properties
    β = medium.ρ*medium.c^2 # medium bulk modulus
    β_eff = real( (1-φ)/β + (φ/η)*sum(s.num_density/(s.ρ*s.c^2) for s in species) )^(-1.0)
    ρ_frac = (φ/η)*sum(s.num_density*(medium.ρ - s.ρ)/(medium.ρ + s.ρ) for s in species)
    ρ_eff = medium.ρ*(1 - ρ_frac)/(1 + ρ_frac)

    return Medium(ρ=ρ_eff, c=sqrt(β_eff/ρ_eff))
end

"calculate effective transmission angle"
function transmission_angle{T<:Number}(k::Complex{T}, k_eff::Complex{T}, θ_inc; tol = 1e-8)
    snell(θ::Array{T}) = norm(k*sin(θ_inc) - k_eff*sin(θ[1] + im*θ[2]))
    result = optimize(snell, [θ_inc,0.]; g_tol= tol^2.0)
    θ_eff = result.minimizer[1] + im*result.minimizer[2]
end

"the reflection from a homogenious halfspace, which is also the low frequency reflection from a particulate material when using the effective_medium."
function reflection_coefficient_halfspace{T<:Number}(incident_medium::Medium{T}, reflect_medium::Medium{T};
        θ_inc::T = zero(T), tol = 1e-6)

    q = real(reflect_medium.c*reflect_medium.ρ/(incident_medium.c*incident_medium.ρ))
    θ_trans = transmission_angle(reflect_medium.c, incident_medium.c,θ_inc; tol=tol)
    R = (q*cos(θ_inc) - cos(θ_trans))/(q*cos(θ_inc) + cos(θ_trans))
end
