"The average reflection coefficients"
function reflection_coefficients(ωs::Union{T,AbstractVector{T}}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number

    k_effs = wavenumber(ωs, medium, species; kws...)
    Rs = [
        reflection_coefficient(ωs[i], EffectiveWave(ωs[i], k_effs[i], medium, species; kws...), medium, species; kws...)
    for i in eachindex(ωs)]

    return Rs
end

"Pairs each ω in ωs with each wave in waves to calculate the refleciton coefficients: reflection_coefficient(ω, wave)"
reflection_coefficients(ωs::AbstractVector{T}, waves::Vector{EffectiveWave{T}}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number =
[ reflection_coefficient(ωs[i], waves[i], medium, species; kws...) for i in eachindex(ωs)]

"The average reflection coefficient"
function reflection_coefficient(ω::T, wave_eff::EffectiveWave{T}, medium::Medium{T}, species::Vector{Specie{T}};
        θin::T = zero(T), x::T = zero(T), kws...) where T<:Number

    k = ω/medium.c
    θ_ref = pi - wave_eff.θ_eff - θin
    S = length(species)
    ho = wave_eff.hankel_order

    kθ = (k*cos(θin) + wave_eff.k_eff*cos(wave_eff.θ_eff))
    R = 2.0im / (k*cos(θin) * kθ)
    R = R*sum(
        exp(im*n*θ_ref + im*x*kθ) * species[l].num_density *
        wave_eff.amplitudes[n+ho+1,l] * Zn(ω,species[l],medium,n)
    for n=-ho:ho, l=1:S)

    return R
end
