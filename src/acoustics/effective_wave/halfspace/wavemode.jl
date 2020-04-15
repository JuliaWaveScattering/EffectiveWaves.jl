
function source_extinction_system(ω::T, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, material::Material{Dim,Halfspace{T,Dim}}) where {T<:Number,Dim}

    θin = transmission_angle(wave_eff, material)
    θ_eff = transmission_angle(psource, material)

    kcos_in = (ω / psource.medium.c) * dot(- conj(material.shape.normal), psource.direction)
    kcos_eff = dot(- conj(material.shape.normal), wave_eff.wavevector)

    ho = wave_eff.basis_order

    extinction_matrix = T(2) .* transpose(vec(
        [
            exp(im*n*(θin - θ_eff)) * number_density(s)
        for n = -ho:ho, s in material.species]
    ))

    forcing = [im * field(psource,zeros(T,Dim),ω) * kcos_in * (kcos_eff - kcos_in)]

    return extinction_matrix, forcing

end

"returns a number a, such that a*As_eff will cancel an incident wave plane wave with incident angle θin."
function scale_mode_amplitudes(ω::T, wave_eff::EffectivePlaneWaveMode{T}, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, material::Material{Dim,Halfspace{T,Dim}}) where {T<:Number,Dim}

    θin = transmission_angle(wave_eff, material)
    θ_eff = transmission_angle(psource, material)

    kcos_in = (ω / psource.medium.c) * dot(- conj(material.shape.normal), psource.direction)
    kcos_eff = dot(- conj(material.shape.normal), wave_eff.wavevector)

    ho = wave_eff.basis_order

    A = T(2) .* transpose(vec(
        [
            exp(im*n*(θin - θ_eff)) * number_density(s)
        for n = -ho:ho, s in material.species]
    ))

    forcing = [im * field(psource,zeros(T,Dim),ω) * kcos_in * (kcos_eff - kcos_in)]

    α = (A * wave_eff.amplitudes[:]) \ forcing

    return α
end

function effective_wavemode(ω::T, k_eff::Complex{T}, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, material::Material{Dim,Halfspace{T,Dim}};
        tol::T = 1e-6,
        extinction_rescale::Bool = true,
        kws...
    )::EffectivePlaneWaveMode{T,Dim} where {T<:AbstractFloat,Dim}

    k = ω/psource.medium.c

    k_vec = transmission_wavevector(k_eff, (ω / psource.medium.c) * psource.direction, material.shape.normal; tol = tol)

    amps = mode_amplitudes(ω, k_eff, psource, material; tol = tol, kws...)
    plane_mode = EffectivePlaneWaveMode(ω, k_vec, amps)
    amps = amps.*scale_mode_amplitudes(ω, plane_mode, psource, material)

    return EffectivePlaneWaveMode(ω, k_vec, amps)
end
