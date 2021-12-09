"""
    reflection_coefficient(ω::T, m_wave::MatchPlaneWaveMode,
        source::PlaneSource, material::Material{2,Halfspace{T,2}};

Calculate the reflection coefficient from a matched wave. This requires using both the discrete part and effective wavemode of the matched wave.
"""
function reflection_coefficient(ω::T, m_wave::MatchPlaneWaveMode{T},
    source::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{2,Halfspace{T,2}};
     # medium::PhysicalMedium, specie::Specie{T};
    kws...) where {T<:AbstractFloat}

    R_eff = reflection_coefficient(ω, m_wave.PlaneWaveModes, source, material; x = m_wave.x_match[end], kws...)
    R_discrete = reflection_coefficient(ω, m_wave.discrete_wave, source, material; kws...)

    return  R_eff + R_discrete
end
