function reflection_coefficient(ω::T, m_wave::MatchWave{T}, medium::PhysicalMedium{T}, specie::Specie{T}; kws...) where {T<:AbstractFloat}

    R_eff = reflection_coefficient(ω, m_wave.effective_waves, medium, [specie]; x = m_wave.x_match[end], kws...)
    R_discrete = reflection_coefficient(ω, m_wave.average_wave, medium, specie; kws...)

    return  R_eff + R_discrete
end
