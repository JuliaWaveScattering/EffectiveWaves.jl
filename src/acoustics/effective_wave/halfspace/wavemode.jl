function boundary_condition_system(ω::T, k_eff::Complex{T}, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, material::Material{Dim,Halfspace{T,Dim}};
        basis_order::Int = 2,
        kws...
    ) where {T<:Number,Dim}

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

    forcing = [im * field(psource,zeros(T,Dim),ω) * kcos_in * (kcos_eff - kcos_in)]

    return extinction_matrix, forcing

end

function effective_wavemode(ω::T, k_eff::Complex{T}, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, material::Material{Dim,Halfspace{T,Dim}};
        tol::T = 1e-6,
        kws...
    )::EffectivePlaneWaveMode{T,Dim} where {T<:AbstractFloat,Dim}

    k = ω/psource.medium.c

    direction = transmission_direction(k_eff, (ω / psource.medium.c) * psource.direction, material.shape.normal; tol = tol)

    amps = eigenvectors(ω, k_eff, psource, material; kws...)

    return EffectivePlaneWaveMode(ω, k_eff, direction, amps)
end
