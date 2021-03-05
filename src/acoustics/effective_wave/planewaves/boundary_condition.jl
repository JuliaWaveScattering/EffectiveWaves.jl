function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigenvectors::Array{Complex{T}}, psource::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{2,Halfspace{T,2}};
        basis_order::Int = 2,
        kws...
    ) where T

    if size(eigenvectors)[end] > 1
        @warn "The effective wavenumber: $k_eff has more than one eigenvector. For plane-waves this case has not been fully implemented"
    end

    scale_number_density = one(T) - one(T) / material.numberofparticles

    direction_eff = transmission_direction(k_eff, (ω / psource.medium.c) * psource.direction, material.shape.normal)

    θ_eff = transmission_angle(direction_eff, material.shape.normal)
    θin  = transmission_angle(psource, material)

    kcos_in = (ω / psource.medium.c) * dot(- conj(material.shape.normal), psource.direction)
    kcos_eff = k_eff * dot(- conj(material.shape.normal), direction_eff)

    extinction_matrix = T(2) .* transpose(vec(
        [
            exp(im*n*(θin - θ_eff)) * scale_number_density * number_density(s)
        for n = -basis_order:basis_order, s in material.species]
    ))

    forcing = [im * field(psource,zeros(T,2),ω) * kcos_in * (kcos_eff - kcos_in)]

    # where extinction_matrix * eigenvectors * α = forcing
    α = (extinction_matrix * reshape(eigenvectors,(:,size(eigenvectors)[end]))) \ forcing

    return α

end

function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigenvectors1::Array{Complex{T}}, eigenvectors2::Array{Complex{T}}, psource::PlaneSource{T,3,1,Acoustic{T,3}}, material::Material{3,Plate{T,3}};
        basis_order::Int = 2,
        kws...
    ) where T

    if size(eigenvectors1)[end] > 1 || size(eigenvectors2)[end] > 1
        @warn "The effective wavenumber: $k_eff has more than one eigenvector. For plane-waves this case has not been fully implemented"
    end

    scale_number_density = one(T) - one(T) / material.numberofparticles

    # First we calculate the innerward pointing normal
    n = material.shape.normal;
    n = n .* sign(real(dot(n,psource.direction)));

    # calculate the distances from the origin to the two faces of the plate
    Z0 = dot(material.shape.normal, material.shape.origin)
    Z1 = Z0 - material.shape.width / T(2)
    Z2 = Z0 + material.shape.width / T(2)

    # next the components of the wavenumbers in the direction of the inward normal
    direction1 = transmission_direction(k_eff, (ω / psource.medium.c) * psource.direction, -n)
    direction2 = transmission_direction(- k_eff, (ω / psource.medium.c) * psource.direction, -n)

    k = (ω / psource.medium.c)
    kz = k * dot(conj(n), psource.direction)
    k1_effz = k_eff * dot(conj(n), direction1)
    k2_effz = k_eff * dot(conj(n), direction2)

    rθφ = cartesian_to_radial_coordinates(psource.direction);

    Ys = spherical_harmonics(basis_order, rθφ[2], rθφ[3]);
    lm_to_n = lm_to_spherical_harmonic_index

    I1 = [ - Ys[lm_to_n(dl,dm)] * im^T(dl+1) * T(2π) * T(1)^dl *
        exp(im * (k_effz - kz)*(Z1 + outer_radius(s2))) *
        scale_number_density * number_density(s2) / (k*kz * (kz - k_effz))
    for dl = 0:basis_order for dm = -dl:dl, s2 in material.species]

    I2 = [ - Ys[lm_to_n(dl,dm)] * im^T(dl+1) * T(2π) * T(1)^dm *
        exp(im * (k_effz - kz)*(Z2 - outer_radius(s2))) *
        scale_number_density * number_density(s2) / (k*kz * (kz + k_effz))
    for dl = 0:basis_order for dm = -dl:dl, s2 in material.species]

    forcing = [ - field(psource,zeros(T,3),ω), T(0)]

    # where extinction_matrix * eigenvectors * α = forcing
    α = (extinction_matrix * reshape(eigenvectors,(:,size(eigvectors)[end]))) \ forcing

    return α

end
