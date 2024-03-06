function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigenvectors::Array{Complex{T}}, psource::PlaneSource{T,2,1,P}, material::Material{Halfspace{T,2}}, ::AbstractPlanarSymmetry{2};
        basis_order::Int = 2,
        kws...
    ) where {T, P <: PhysicalMedium{2,1}}

    if size(eigenvectors)[end] > 1
        @warn "The effective wavenumber: $k_eff has more than one eigenvector. For plane-waves this case has not been fully implemented"
    end

    if psource.medium != material.microstructure.medium @error mismatched_medium end

    scale_number_density = one(T)

    direction_eff = transmission_direction(k_eff, (ω / psource.medium.c) * psource.direction, material.shape.normal)

    θ_eff = transmission_angle(direction_eff, material.shape.normal)
    θin  = transmission_angle(psource, material)

    kcos_in = (ω / psource.medium.c) * dot(- conj(material.shape.normal), psource.direction)
    kcos_eff = k_eff * dot(- conj(material.shape.normal), direction_eff)

    extinction_matrix = T(2) .* transpose(vec(
        [
            exp(im*n*(θin - θ_eff)) * scale_number_density * number_density(s)
        for n = -basis_order:basis_order, s in material.microstructure.species]
    ))

    forcing = [im * field(psource,zeros(T,2),ω) * kcos_in * (kcos_eff - kcos_in)]

    # where extinction_matrix * eigenvectors * α = forcing
    α = (extinction_matrix * reshape(eigenvectors,(:,size(eigenvectors)[end]))) \ forcing

    return α

end

function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, psource::PlaneSource{T,3,1,P}, material::Material{Halfspace{T,3}}, ::AbstractPlanarSymmetry{3};
        basis_order::Int = 2,
        kws...
    ) where {T, P <: PhysicalMedium{3,1}}

    if size(eigvectors)[end] > 1
        @warn "The effective wavenumber: $k_eff has more than one eigenvector. For plane-waves this case has not been fully implemented"
    end
    if psource.medium == material.microstructure.medium  

        # First we calculate the outward pointing normal
        n = material.shape.normal;
        n = - n .* sign(real(dot(n,psource.direction)));

        # calculate the distances from the origin to the face of the halfspace
        Z0 = dot(-n, material.shape.origin)

        # next the components of the wavenumbers in the direction of the inward normal
        direction = transmission_direction(k_eff, (ω / psource.medium.c) * psource.direction, n)

        k = (ω / psource.medium.c)
        kz = k * dot(-conj(n), psource.direction)
        k_effz = k_eff * dot(-conj(n), direction)
        rθφ = cartesian_to_radial_coordinates(psource.direction);

        Ys = spherical_harmonics(basis_order, rθφ[2], rθφ[3]);  
        lm_to_n = lm_to_spherical_harmonic_index

        I1(F::Array{Complex{T}}, k_effz::Complex{T}) = sum(
            - F[lm_to_n(dl,dm),j,p] * Ys[lm_to_n(dl,dm)] * im^T(dl+1) * T(2π) * T(-1)^dl *
            exp(im * (k_effz - kz)*(Z0 + outer_radius(material.microstructure.species[j]))) * number_density(material.microstructure.species[j]) / (k*kz * (kz - k_effz))
        for p = 1:size(F,3), dl = 0:basis_order for dm = -dl:dl, j in eachindex(material.microstructure.species))

        forcing = - field(psource,zeros(T,3),ω)

        α = I1(eigvectors,k_effz) \ forcing

        return α

    else
        
        # Setting parameters
        c = psource.medium.c
        c0 = material.microstructure.medium.c
        k = ω / c
        k0 = ω / c0
        species = material.microstructure.species
        rs = outer_radius.(species)
        nf = number_density(material.microstructure.species)
        Fs = eigvectors
        M = size(Fs,3)

        if (c0 != real(c0) || c != real(c))
            @warn "The case of complex wavespeed has not been fully implemented for two medium case, using real part of wavespeeds instead."
        end

        # First we calculate the outward pointing normal
        n = material.shape.normal;
        n = - n .* sign(real(dot(n,psource.direction)));

        # calculate the distances from the origin to the face of the halfspace
        Z0 = dot(-n, material.shape.origin)

        # next the components of the wavenumbers in the direction of the inward normal
        kz = k * dot(-conj(n), psource.direction)
        k0z = sqrt(k0^2 - (k^2 - kz^2))

        # Needs adjustments for complex wavespeeds
        direction0 = real(k) * (psource.direction - dot(-conj(n), psource.direction) * psource.direction)
        direction0 += real(k0z) * dot(-conj(n), psource.direction) * psource.direction
        direction0 /= norm(direction0)

        direction = transmission_direction(k_eff, k * psource.direction, n)
        
        k_effz = k_eff * dot(-conj(n), direction)

        γ0 = kz/k0z
        ζR = (kz - k0z) / (kz + k0z)
        ζT = 2k0z / (kz + k0z)
        
        rθφ = cartesian_to_radial_coordinates(direction0)
        # spherical_harmonics() needs adjustments for complex wavespeeds
        Ys = spherical_harmonics(basis_order, rθφ[2], rθφ[3])
        lm_to_n = lm_to_spherical_harmonic_index

        particle_contribution = sum(
            Fs[lm_to_n(dl, dm), j, p] * nf[j] * (-1im)^dl * Ys[lm_to_n(dl, dm)] * exp(1im * (k_effz - k0z) * (Z0 + rs[j]))
        for p = 1:M, dl = 0:basis_order for dm = -dl:dl, j in eachindex(species))

        particle_contribution *= (2pi * 1im) * (k_effz + k0z) / (k0 * k0z * (k0^2 - k_eff^2))
        
        rθφ_n = cartesian_to_radial_coordinates(-conj(n))
        Ys_n = spherical_harmonics(basis_order, rθφ_n[2], rθφ_n[3])

        wall_contribution = sum(
            -Fs[lm_to_n(dl,dm),j,p] * nf[j] * 1im^dl * Ys_n[lm_to_n(dl, dm)] * exp(1im * (k_effz + k0z) * (Z0 + rs[j]))
        for p = 1:M, dl = 0:basis_order for dm = -dl:dl, j in eachindex(species))

        wall_contribution *= ζR * 2 * pi / (1im * k0 * k0z * (k_effz + k0z))

        M_dot = particle_contribution + wall_contribution
        
        forcing = γ0 * ζT * field(psource, zeros(T, 3), ω)

        α = M_dot \ forcing

        return α
    end

end

function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors1::Array{Complex{T}}, eigvectors2::Array{Complex{T}}, psource::PlaneSource{T,3,1,P}, material::Material{Plate{T,3}}, ::AbstractPlanarSymmetry{3};
        basis_order::Int = 2,
        kws...
    ) where {T, P <: PhysicalMedium{3,1}}

    if psource.medium != material.microstructure.medium @error mismatched_medium end

    if size(eigvectors1)[end] > 1 || size(eigvectors2)[end] > 1
        @warn "The effective wavenumber: $k_eff has more than one eigenvector. For plane-waves this case has not been fully implemented"
    end

    scale_number_density = one(T)

    # First we calculate the outward pointing normal
    n = material.shape.normal;
    n = - n .* sign(real(dot(n,psource.direction)));

    # calculate the distances from the origin to the two faces of the plate
    Z0 = dot(-n, material.shape.origin)
    Z1 = Z0 - material.shape.width / 2
    Z2 = Z0 + material.shape.width / 2

    # next the components of the wavenumbers in the direction of the inward normal
    k = (ω / psource.medium.c)
    k1 = k_eff;
    k2 = -k_eff;

    direction1 = transmission_direction(k1, k * psource.direction, n)
    direction2 = direction1

    kz = k * dot(-conj(n), psource.direction)
    k1_effz = k1 * dot(-conj(n), direction1)
    k2_effz = k2 * dot(-conj(n), direction2)

    rθφ = cartesian_to_radial_coordinates(psource.direction);

    Ys = spherical_harmonics(basis_order, rθφ[2], rθφ[3]);
    lm_to_n = lm_to_spherical_harmonic_index

    I1(F::Array{Complex{T}}, k_effz::Complex{T}) = sum(
        - F[lm_to_n(dl,dm),j,p] * Ys[lm_to_n(dl,dm)] * im^T(dl+1) * T(2π) * T(-1)^dl *
        exp(im * (k_effz - kz)*(Z1 + outer_radius(material.microstructure.species[j]))) *
        scale_number_density * number_density(material.microstructure.species[j]) / (k*kz * (kz - k_effz))
    for p = 1:size(F,3), dl = 0:basis_order for dm = -dl:dl, j in eachindex(material.microstructure.species))

    I2(F::Array{Complex{T}}, k_effz::Complex{T}) = sum(
        - F[lm_to_n(dl,dm),j,p] * Ys[lm_to_n(dl,dm)] * im^T(dl+1) * T(2π) * T(-1)^dm *
        exp(im * (k_effz + kz)*(Z2 - outer_radius(material.microstructure.species[j]))) *
        scale_number_density * number_density(material.microstructure.species[j]) / (k*kz * (kz + k_effz))
    for p = 1:size(F,3), dl = 0:basis_order for dm = -dl:dl, j in eachindex(material.microstructure.species))

    MIs = [I1(eigvectors1,k1_effz) I1(eigvectors2,k2_effz);
           I2(eigvectors1,k1_effz) I2(eigvectors2,k2_effz)]

    forcing = [ - field(psource,zeros(T,3),ω), T(0)]

    α = MIs \ forcing

    return α

end
