"The effective low frequency bulk modulus and density of a material filled with particles"
function effective_medium(medium::Acoustic{T,Dim}, species::Species{T,Dim};
        numberofparticles::Number = Inf
    ) where {T<:AbstractFloat,Dim}

    scale_number_density = one(T) - one(T) / numberofparticles
    φ = scale_number_density * sum(volume_fraction.(species))

    # calculate effective properties
    β = medium.ρ*medium.c^2 # medium bulk modulus
    if abs(β) == zero(T) || abs(prod(s -> s.particle.medium.ρ * s.particle.medium.c^2, species)) == zero(T)
        β_eff = zero(T)
    else
        β_eff = ((one(T)-φ) / β + scale_number_density * sum(s.volume_fraction / (s.particle.medium.ρ * s.particle.medium.c^2) for s in species))^(-one(T))
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

    return Acoustic(Dim; ρ = ρ_eff, c = sqrt(β_eff/ρ_eff))
end

"
reflection_coefficient(PlaneSource, Acoustic[, Halfspace = Halfspace(-psource.direction)])

caculates the reflection coefficient from a homogenious halfspace (assumed to direct incidence if not given), which is also the low frequency reflection from a particulate material when using the effective_medium."
function reflection_coefficient(ω::T, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, reflect_medium::Acoustic{T,Dim}, halfspace::Halfspace{T,Dim} = Halfspace(-psource.direction)) where {T<:AbstractFloat,Dim}

    k_in = ω / psource.medium.c
    k_r = if abs(reflect_medium.c) == zero(T)
        T(Inf) + zero(T) * im
    else
        ω / reflect_medium.c
    end

    θin = transmission_angle(psource, halfspace)

    vtran = transmission_direction(k_r, k_in * psource.direction, halfspace.normal)
    θtran = transmission_angle(vtran, halfspace.normal)

    q = real(reflect_medium.c*reflect_medium.ρ / (psource.medium.c*psource.medium.ρ))
    f = field(psource)(halfspace.origin,ω)

    R = (f^2 / psource.amplitude[1]) * (q*cos(θin) - cos(θtran)) / (q*cos(θin) + cos(θtran))
end

"""
    planewave_coefficients(psource::PlaneSource, reflect_medium::Acoustic, plate::Plate)

Calculates the coefficients, or amplityudes, of the plane waves scattering from and transmitted inside the plate, for the plane wave source `psource` and the `plate`. The function returns `[R, T, P1, P2]` where `R` is the reflection coefficient, `T` is the coefficient of the transmitted wave, and `P1` (`P2`) are the amplitudes of the wave travling forward (backward) inside the plate.
"""
function planewave_coefficients(ω::T, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, plate_medium::Acoustic{T,Dim}, plate::Plate{T,Dim}) where {T<:AbstractFloat,Dim}

    k0 = ω / psource.medium.c
    k1 = if abs(plate_medium.c) == zero(T)
        T(Inf) + zero(T) * im
    else
        ω / plate_medium.c
    end

    # make the normal outward pointing
    n = plate.normal / norm(plate.normal)
    if real(dot(n,psource.direction)) > 0
        n = -n
    end

    θ0 = transmission_angle(psource, plate)

    v1 = transmission_direction(k1, k0 .* psource.direction, n)
    θ1 = transmission_angle(v1, n)

    C0 = k0 * cos(θ0) / psource.medium.ρ
    C1 = k1 * cos(θ1) / plate_medium.ρ

    Z0 = dot(-n,plate.origin - psource.position)
    Z1 = Z0 - plate.width / T(2)
    Z2 = Z0 + plate.width / T(2)

    UR = (C0 - C1)*(C0 + C1) * exp(2im*k0*Z1*cos(θ0)) * (exp(2im*k1*Z1*cos(θ1)) - exp(2im*k1*Z2*cos(θ1))) /
        ((C0 + C1)^2 * exp(2im*k1*Z1*cos(θ1)) - (C0 - C1)^2 * exp(2im*k1*Z2*cos(θ1)))

    UP1 = 2C0*(C0 + C1) * exp(im*Z1*(k0*cos(θ0) + k1*cos(θ1))) /
        ((C0 + C1)^2*exp(2im*k1*Z1*cos(θ1)) - (C0 - C1)^2*exp(2im*k1*Z2*cos(θ1)))

    UP2 = 2*C0*(C0 - C1)*exp(im*(k0*Z1*cos(θ0) + k1*(Z1 + 2*Z2)*cos(θ1))) /
        (-((C0 + C1)^2*exp((2*im)*k1*Z1*cos(θ1))) + (C0 - C1)^2*exp(2im*k1*Z2*cos(θ1)))

    UT = 4*C0*C1*exp(im*(k0*(Z1 - Z2)*cos(θ0) + k1*(Z1 + Z2)*cos(θ1))) /
    ((C0 + C1)^2*exp(2im*k1*Z1*cos(θ1)) - (C0 - C1)^2*exp(2im*k1*Z2*cos(θ1)))

    return [UR, UT, UP1, UP2]
end
