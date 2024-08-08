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
        ρ = psource.medium.ρ
        ρ0 = material.microstructure.medium.ρ
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
        direction_k0 = real(k) * (psource.direction - dot(-conj(n), psource.direction) * psource.direction)
        direction_k0 = direction_k0 + real(k0z) * dot(-conj(n), psource.direction) * psource.direction
        direction_k0 = direction_k0/norm(direction_k0)

        direction = transmission_direction(k_eff, k * psource.direction, n)
        
        k_effz = k_eff * dot(-conj(n), direction)

        γ0 = (ρ0 * kz) /(ρ * k0z)
        ζR = (ρ0 * kz - ρ * k0z) / (ρ0 * kz + ρ * k0z)
        ζT = 2ρ * k0z / (ρ0 * kz + ρ * k0z)
        
        rθφ = cartesian_to_radial_coordinates(direction_k0)
        # spherical_harmonics() needs adjustments for complex wavespeeds
        Ys = spherical_harmonics(basis_order, rθφ[2], rθφ[3])
        lm_to_n = lm_to_spherical_harmonic_index

        particle_contribution = sum(
            Fs[lm_to_n(dl, dm), j, p] * nf[j] * (-1im)^dl * Ys[lm_to_n(dl, dm)] * exp(1im * (k_effz - k0z) * (Z0 + rs[j]))
        for p = 1:M, dl = 0:basis_order for dm = -dl:dl, j in eachindex(species))

        particle_contribution *= (2pi * 1im) * (k_effz + k0z) / (k0 * k0z * (k0^2 - k_eff^2))

        wall_contribution = sum(
            -Fs[lm_to_n(dl,dm),j,p] * nf[j] * 1im^dl * Ys[lm_to_n(dl, dm)] * exp(1im * (k_effz + k0z) * (Z0 + rs[j]))
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

        if size(eigvectors1)[end] > 1 || size(eigvectors2)[end] > 1
        @warn "The effective wavenumber: $k_eff has more than one eigenvector. For plane-waves this case has not been fully implemented"
    end

    if psource.medium == material.microstructure.medium

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
    
    else

        # Setting parameters
        ρ = psource.medium.ρ
        ρ0 = material.microstructure.medium.ρ
        c = psource.medium.c
        c0 = material.microstructure.medium.c
        k = ω / c
        k0 = ω / c0
        species = material.microstructure.species
        rs = outer_radius.(species)
        nf = number_density.(species)
        F_p = eigvectors1
        F_m = eigvectors2
        Z = material.shape.width
        G = field(psource, zeros(T, 3), ω)

        if (c0 != real(c0) || c != real(c))
            @warn "The case of complex wavespeed has not been fully implemented for two medium case, using real part of wavespeeds instead."
        end

        # First we calculate the outward pointing normal
        n = material.shape.normal
        if norm(dot(n, psource.direction)) != norm(n) * norm(psource.direction)
            @warn "Only normal incidence has been implemented yet for a two medium plate case. Only normal direction of the plane wave will contribute to reflection and transmission."
        end
        n = -n .* sign(real(dot(n, psource.direction)))

        # calculate the distances from the origin to the two faces of the plate
        Z0 = dot(-n, material.shape.origin)
        Z1 = Z0 - Z / 2
        Z2 = Z0 + Z / 2

        # next the components of the wavenumbers in the direction of the inward normal
        kz = k * dot(-conj(n), psource.direction)
        k0z = sqrt(k0^2 - (k^2 - kz^2))

        # Needs adjustments for complex wavespeeds
        direction_k0 = real(k) * (psource.direction - dot(-conj(n), psource.direction) * psource.direction)
        direction_k0 = direction_k0 + real(k0z) * dot(-conj(n), psource.direction) * psource.direction
        direction_k0 = direction_k0 / norm(direction_k0)

        direction = transmission_direction(k_eff, k * psource.direction, n)

        k_effz_p = k_eff * dot(-conj(n), direction)
        k_effz_m = -k_effz_p

        # Needed coefficients
        Δ = 2 * k * k0 * ρ * ρ0 * cos(k0*Z) - 1im * ((k0 * ρ)^2 + (k * ρ0)^2) * sin(k0 * Z)
        D_p = k0 * ρ * (k0 * ρ + k * ρ0) / Δ
        D_m = k0 * ρ * (k0 * ρ - k * ρ0) / Δ
        D_1 = ((k0 * ρ)^2 - (k * ρ0)^2) / 2Δ
        D_2 = (k0 * ρ - k * ρ0)^2 / 2Δ
        γ0 = k * ρ0 / (k0 * ρ)

        rθφ = cartesian_to_radial_coordinates(direction)

        Ys = spherical_harmonics(basis_order, rθφ[2], rθφ[3])
        lm_to_n = lm_to_spherical_harmonic_index

        Pr(x::Complex{T}, y::Complex{T}, r::T) = exp(1im * (x + y) * Z / 2) * sin((x + y) * (Z / 2 - r)) / (x + y)

        Bp(F::Array{Complex{T}}, p_or_m::Int) = (4π / (k0^2)) * sum(
            1im^(-T(l)) * Ys[lm_to_n(l, m)] * Pr(p_or_m * k_eff, -k0, rs[j]) * F[lm_to_n(l, m), j, p] * nf[j]
        for p = 1:size(F, 3), l = 0:basis_order for m = -l:l, j in eachindex(species))

        Bm(F::Array{Complex{T}}, p_or_m::Int) = (4π / (k0^2)) * sum(
            1im^(T(l)) * Ys[lm_to_n(l, m)] * Pr(p_or_m * k_eff, k0, rs[j]) * F[lm_to_n(l, m), j, p] * nf[j]
        for p = 1:size(F, 3), l = 0:basis_order for m = -l:l, j in eachindex(species))

        Iu(F::Array{Complex{T}}, p_or_m::Int) = (2π / (k0^2)) * sum(
            1im^(-T(dl)) * Ys[lm_to_n(dl, dm)] * exp(1im * (p_or_m * k_eff - k0) * rs[j]) * F[lm_to_n(dl, dm), j, p] * nf[j] / (k_eff - p_or_m * k0)
        for p = 1:size(F, 3), dl = 0:basis_order for dm = -dl:dl, j in eachindex(species))

        Il(F::Array{Complex{T}}, p_or_m::Int) = (2π / (k0^2)) * sum(
            1im^(T(dl)) * Ys[lm_to_n(dl, dm)] * exp(1im * (p_or_m * k_eff + k0) * (Z - rs[j])) * F[lm_to_n(dl, dm), j, p] * nf[j] / (k_eff + p_or_m * k0)
        for p = 1:size(F, 3), dl = 0:basis_order for dm = -dl:dl, j in eachindex(species))

        MI11 = D_2 * Bp(F_p, 1) * exp(1im * k0 * Z) + D_1 * Bm(F_p, 1) * exp(-1im * k0 * Z) + 1im * Iu(F_p, 1)
        MI12 = D_2 * Bp(F_m, -1) * exp(1im * k0 * Z) + D_1 * Bm(F_m, -1) * exp(-1im * k0 * Z) + 1im * Iu(F_m, -1)
        MI21 = (D_1 * Bp(F_p, 1) + D_2 * Bm(F_p, 1)) * exp(1im * k0 * Z) - 1im * Il(F_p, 1)
        MI22 = (D_1 * Bp(F_m, -1) + D_2 * Bm(F_m, -1)) * exp(1im * k0 * Z) - 1im * Il(F_m, -1)

        MIs = [MI11 MI12;
               MI21 MI22]

        forcing = -G * γ0 * [D_p * exp(-1im * k0 * Z), D_m * exp(1im * k0 * Z)]

        α = MIs \ forcing

        return α

    end    

end
