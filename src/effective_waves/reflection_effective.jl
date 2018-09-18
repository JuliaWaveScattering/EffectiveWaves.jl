"The average reflection coefficients"
function reflection_coefficients(ωs::Union{T,AbstractVector{T}}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number

    Rs = map(ωs) do ω
        k_effs = wavenumbers(ω, medium, species; kws...)
        w_effs = [EffectiveWave(ω, k_eff, medium, species; kws...) for k_eff in k_effs]
        reflection_coefficient(ω, w_effs, medium, species; kws...)
    end

    return Rs
end

"Calculates the reflection coefficient from each wave in waves and then sums the results."
reflection_coefficient(ω::T, waves::Vector{EffectiveWave{T}}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number =
    sum(w -> reflection_coefficient(ω, w, medium, species; kws...), waves)

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
        wave_eff.amplitudes[n+ho+1,l]
    for n=-ho:ho, l=1:S)

    return R
end
