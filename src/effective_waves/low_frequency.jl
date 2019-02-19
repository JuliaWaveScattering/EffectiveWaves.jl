"the effective low frequency bulk modulus and density of a material filled with particles"
function effective_medium(medium::Medium{T}, species::Array{Specie{T}}) where T<:Number
    φ = sum(volume_fraction.(species))

    # calculate effective properties
    β = medium.ρ*medium.c^2 # medium bulk modulus
    if abs(β) == zero(T) || abs(prod(s -> s.ρ*s.c^2, species)) == zero(T)
        β_eff = zero(T)
    else
        β_eff = one(T)/((1-φ)/β + sum(volume_fraction(s)/(s.ρ*s.c^2) for s in species))
    end
    ρ_frac = sum(volume_fraction(s)*(medium.ρ - s.ρ)/(medium.ρ + s.ρ) for s in species)
    ρ_eff = medium.ρ*(1 - ρ_frac)/(1 + ρ_frac)

    return Medium(ρ=ρ_eff, c=sqrt(β_eff/ρ_eff))
end

"the reflection from a homogenious halfspace, which is also the low frequency reflection from a particulate material when using the effective_medium."
function reflection_coefficient_halfspace(incident_medium::Medium{T}, reflect_medium::Medium{T};
        θin::T = zero(T), tol = 1e-6) where T<:Number

    q = real(reflect_medium.c*reflect_medium.ρ/(incident_medium.c*incident_medium.ρ))
    θ_trans = transmission_angle(reflect_medium.c, incident_medium.c, θin; tol=tol)
    R = (q*cos(θin) - cos(θ_trans))/(q*cos(θin) + cos(θ_trans))
end
