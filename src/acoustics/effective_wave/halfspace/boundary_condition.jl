function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigenvectors::Array{Complex{T}}, psource::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{2,Halfspace{T,2}};
        basis_order::Int = 2,
        kws...
    ) where T

    if size(eigenvectors,2) > 1
        @warn "The effective wavenumber: $k_eff has more than one eigenvector. For plane-waves this case has not been fully implemented"
    end

    direction = transmission_direction(k_eff, (ω / psource.medium.c) * psource.direction, material.shape.normal)

    θin = transmission_angle(direction, material.shape.normal)
    θ_eff = transmission_angle(psource, material)

    kcos_in = (ω / psource.medium.c) * dot(- conj(material.shape.normal), psource.direction)
    kcos_eff = k_eff * dot(- conj(material.shape.normal), direction)

    extinction_matrix = T(2) .* transpose(vec(
        [
            exp(im*n*(θin - θ_eff)) * number_density(s)
        for n = -basis_order:basis_order, s in material.species]
    ))

    forcing = [im * field(psource,zeros(T,2),ω) * kcos_in * (kcos_eff - kcos_in)]

    # where extinction_matrix * eigenvectors * α = forcing
    α = (extinction_matrix * eigenvectors) \ forcing

    return α

end
