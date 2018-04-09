function effective_material_properties{T<:Number}(medium::Medium{T}, species::Array{Specie{T}};
        tol = 1e-6, θ_inc::T = 0.0,
        kws...)
    # total number density
    η = sum(s.num_density for s in species)
    # total volume fraction
    φ = sum(s.num_density*(s.r^2)*pi for s in species)
    β = medium.ρ*medium.c^2
    β_eff = real( (1-φ)/β + (φ/η)*sum(s.num_density/(s.ρ*s.c^2) for s in species) )^(-1.0)
    ρ_frac = (φ/η)*sum(s.num_density*(medium.ρ - s.ρ)/(medium.ρ + s.ρ) for s in species)
    ρ_eff = medium.ρ*(1 - ρ_frac)/(1 + ρ_frac)
    return (β_eff, ρ_eff)
end

"the reflection from a homogenious halfspace, which is also the low frequency reflection from a particulate material when using the effective_material_properties."
function reflection_coefficient_halfspace{T<:Number}(incident_medium::Medium{T}, θ_inc::T, reflect_medium::Medium{T})
    q = reflect_medium.c*reflect_medium.ρ/(incident_medium.c*incident_medium.ρ)
    function snell!(F,x)
        F[1] = real(reflect_medium.c*sin(θ_inc) - incident_medium.c*sin(x[1]))
    end
end
