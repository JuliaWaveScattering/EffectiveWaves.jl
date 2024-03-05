effective_medium(medium::PhysicalMedium, sps::Species; kws...) = effective_medium(Microstructure(medium, sps); kws...)

"The effective low frequency bulk modulus and density of a material filled with particles"
function effective_medium(micro::Microstructure{Dim};
        # numberofparticles::Number = Inf
    ) where Dim

    medium = micro.medium
    T = typeof(medium.ρ)

    species = micro.species
    scale_number_density = one(T)
    φ = scale_number_density * volume_fraction(species)

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

"""
reflection_coefficient(PlaneSource, Acoustic[, Halfspace = Halfspace(-psource.direction)])

calculates the reflection coefficient from a homogenious halfspace (assumed to direct incidence if not given), which is also the low frequency reflection from a particulate material when using the effective_medium.
"""
reflection_coefficient(ω::T, source::PlaneSource{T}, reflect_medium::Acoustic{T}, halfspace::Halfspace{T} = Halfspace(-source.direction)) where T<:AbstractFloat = reflection_transmission_coefficients(ω, source, reflect_medium, halfspace)[1]

"""
reflection_transmission_coefficients(PlaneSource, Acoustic[, Halfspace = Halfspace(-psource.direction)])

Calculates the reflection and transmission coefficient, R and T, from a homogenious halfspace, which is also the low frequency reflection from a particulate material when using the effective_medium. 

Let ``\\mathbf k_R`` and ``\\mathbf k_T`` be the wave vectors of the refected and transmitted waves ``u_R(x)`` and ``u_T(x)``, then can describe these waves as

``u_R(x) = R \\mathrm e^{i \\mathbf k_R \\cdot (\\mathbf x - \\mathbf x_1)}``
``u_T(x) = T \\mathrm e^{i \\mathbf k_T \\cdot (\\mathbf x - \\mathbf x_1)}``
"""
function reflection_transmission_coefficients(ω::T, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, reflect_medium::Acoustic{T,Dim}, halfspace::Halfspace{T,Dim} = Halfspace(-psource.direction)) where {T<:AbstractFloat,Dim}

    k_in = ω / psource.medium.c
    k_r = if abs(reflect_medium.c) == zero(T)
        T(Inf) + zero(T) * im
    else
        ω / reflect_medium.c
    end

    # make sure the normal outward pointing
    n = halfspace.normal / norm(halfspace.normal)
    if real(dot(n,psource.direction)) > 0
        n = -n
    end

    θin = transmission_angle(psource.direction, n)

    vtran = transmission_direction(k_r, k_in * psource.direction, n)
    θtran = transmission_angle(vtran, n)

    q = reflect_medium.c*reflect_medium.ρ / (psource.medium.c*psource.medium.ρ)

    Z = dot(n,halfspace.origin)

    R1 = (exp(2im * k_in * Z * cos(θin)) * field(psource)(zeros(T,Dim),ω)) * 
        (q*cos(θin) - cos(θtran)) / (q*cos(θin) + cos(θtran))
    R1 = R1 * exp(im * k_in * Z * cos(θin)) # assuming a reflected wave  = R * exp(i kR.(X-X1)) where kR.X1 =  k_in * Z * cos(θin)

    T1 = exp(im * Z * (k_in * cos(θin) - k_r * cos(θtran))) * field(psource)(zeros(T,Dim),ω) * 2*q*cos(θin) / (q*cos(θin) + cos(θtran))
    T1 = T1 * exp(im * k_r * Z * cos(θtran)) # assuming a transmitter wave  = T * exp(i kT.(X-X1)) where kT.X1 =  k_r * Z * cos(θtran)

    return [R1,T1]
end

"""
    planewave_amplitudes(psource::PlaneSource, reflect_medium::Acoustic, plate::Plate)

Calculates the coefficients, or amplitudes, of the plane waves scattering from and transmitted inside the plate, for the plane wave source `psource` and the `plate`. The function returns `[R, T, P1, P2]` where `R` is the reflection coefficient, `T` is the coefficient of the transmitted wave, and `P1` (`P2`) are the amplitudes of the wave travling forward (backward) inside the plate. 

In more detail, let ``\\mathbf k_R, \\mathbf k_{P_1}, \\mathbf k_{P_2},`` and ``\\mathbf k_T`` be the wave vectors of the refected wave ``u_R(x)``, the wave transmitted into the plate ``u_{P_1}(x)``, the wave in the plate ``u_{P_2}(x)`` (which is the result of ``u_{P_1}(x)`` being reflected from the second plate face) and the wave transmitted to the other side of the plate ``u_T(x)``. To simplify, let us consider a plate occupying the region ``z_1<z<z_2`` and an incident wave vector ``\\mathbf k = [k_x,0,k_z]``, although we note the code works for any plate orientation and source. With these choices we have that ``\\mathbf k_R = [k_x,0,-k_z]``, ``\\mathbf k_{P_1} = [``\\kappa_x,0,\\kappa_z]`` and ``\\mathbf k_{P_2} = [``\\kappa_x,0,-\\kappa_z]``, where ``\\kappa_x`` and ``\\kappa_z`` can be deduced from Snells law (i.e. the boundary conditions). The waves involved are then given by

UR * exp(im*KR.(X-X1))

``u_R(x) = R \\mathrm e^{i \\mathbf k_R \\cdot (\\mathbf x - \\mathbf x_1)}``
``u_{P_1}(x) = P_1 \\mathrm e^{i \\mathbf k_{P_1} \\cdot (\\mathbf x - \\mathbf x_1)}``
``u_{P_2}(x) = P_2 \\mathrm e^{i \\mathbf k_{P_2} \\cdot (\\mathbf x - \\mathbf x_2)}``
``u_T(x) = T \\mathrm e^{i \\mathbf k_T \\cdot (\\mathbf x - \\mathbf x_2)}``

where ``\\mathbf x_1 = [0,0,z_1]`` and ``\\mathbf x_2 = [0,0,z_2]`` for the simplifying choices made above. We note that the phase factors chosen for each of these waves is important and not arbitrary: when the waves attenuate in the direction of propagation they also grow in the opposite direction. The phase factors are chosen so that these waves have reasonable values near the boundaries.
"""
function planewave_amplitudes(ω::T, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, plate_medium::Acoustic{T,Dim}, plate::Plate{T,Dim}) where {T<:AbstractFloat,Dim}

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

    Z0 = dot(-n,plate.origin)
    Z1 = Z0 - plate.width / T(2)
    Z2 = Z0 + plate.width / T(2)

    denom = ((C0 + C1)^2 * exp(2im*k1*Z1*cos(θ1)) - (C0 - C1)^2 * exp(2im*k1*Z2*cos(θ1)))

    UR = (C0 - C1)*(C0 + C1) * exp(2im*k0*Z1*cos(θ0)) * (exp(2im*k1*Z1*cos(θ1)) - exp(2im*k1*Z2*cos(θ1))) / denom
    UR = UR * exp(im*k0*Z1*cos(θ0)) # we assume the reflected wave = UR * exp(im*KR.(X-X1)), where KR.X1 = im*k0*Z1*cos(θ0) = k0 .* [sin(θ0),0,-cos(θ0)], and X1 = [0.0,0.0,Z1], so that KR.X1 = im*k0*Z1*cos(θ0).

    UT = 4*C0*C1*exp(im*(k0*(Z1 - Z2)*cos(θ0) + k1*(Z1 + Z2)*cos(θ1))) / denom
    UT = UT * exp(im*k0*Z2*cos(θ0))

    UP1 = 2C0*(C0 + C1) * exp(im*Z1*(k0*cos(θ0) + k1*cos(θ1))) / denom
    UP1 = UP1 * exp(im*k1*Z1*cos(θ1))

    UP2 = - 2*C0*(C0 - C1)*exp(im*(k0*Z1*cos(θ0) + k1*(Z1 + 2*Z2)*cos(θ1))) / denom
    UP2 = UP2 * exp(im*k1*Z2*cos(θ1))

    return [UR, UT, UP1, UP2] * field(psource)(zeros(T,3),ω)
end
