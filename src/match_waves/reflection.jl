function reflection_coefficient(ω::T, m_wave::MatchPlaneWaveMode{T},
    source::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{2,Halfspace{T,2}};
     # medium::PhysicalMedium{T}, specie::Specie{T};
    kws...) where {T<:AbstractFloat}

    R_eff = reflection_coefficient(ω, m_wave.effective_wavemodes, source, material; x = m_wave.x_match[end], kws...)
    R_discrete = reflection_coefficient(ω, m_wave.discrete_wave, source, material; kws...)

    return  R_eff + R_discrete
end
