"The effective low frequency bulk modulus and density of a material filled with particles"
function effective_medium(medium::Acoustic{T,Dim}, species::Species{T,Dim},
        numberofparticles::Number = Inf
    ) where {T<:AbstractFloat,Dim}
    φ = sum(EffectiveWaves.volume_fraction.(species))

    scale_number_density = one(T)
    if numberofparticles > 1
        scale_number_density = one(T) - one(T) / numberofparticles
    end
    φ = scale_number_density * sum(volume_fraction.(species))


    # calculate effective properties
    β = medium.ρ*medium.c^2 # medium bulk modulus
    if abs(β) == zero(T) || abs(prod(s -> s.particle.medium.ρ * s.particle.medium.c^2, species)) == zero(T)
        β_eff = zero(T)
    else
        β_eff = one(T) / ((1-φ) / β + scale_number_density * sum(s.volume_fraction / (s.particle.medium.ρ * s.particle.medium.c^2) for s in species))
    end

    ρ_frac = scale_number_density * sum(
        map(species) do s
            if s.particle.medium.ρ == Inf
                - s.volume_fraction / T(Dim - 1)
            else
                s.volume_fraction * (medium.ρ - s.particle.medium.ρ)/(medium.ρ + T(Dim - 1) * s.particle.medium.ρ)
            end
        end
    )

    ρ_eff = medium.ρ*(1 - ρ_frac)/(1 + T(Dim - 1) * ρ_frac)

    return Acoustic(Dim; ρ=ρ_eff, c=sqrt(β_eff/ρ_eff))
end

"
reflection_coefficient(PlaneSource, Acoustic[, Halfspace = Halfspace(-psource.direction)])

caculates the reflection coefficient from a homogenious halfspace (assumed to direct incidence if not given), which is also the low frequency reflection from a particulate material when using the effective_medium."
function reflection_coefficient(psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, reflect_medium::Acoustic{T,Dim}, halfspace::Halfspace{T,Dim} = Halfspace(-psource.direction)) where {T<:AbstractFloat,Dim}

    k_in = T(1) / psource.medium.c
    k_r = if abs(reflect_medium.c) == zero(T)
        T(Inf) + zero(T) * im
    else
        T(1)/reflect_medium.c
    end

    θin = transmission_angle(psource, halfspace)

    vtran = transmission_direction(k_r, k_in * psource.direction, halfspace.normal)
    θtran = transmission_angle(vtran, halfspace.normal)

    q = real(reflect_medium.c*reflect_medium.ρ/(psource.medium.c*psource.medium.ρ))
    R = (q*cos(θin) - cos(θtran))/(q*cos(θin) + cos(θtran))
end
