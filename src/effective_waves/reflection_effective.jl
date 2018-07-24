"The average reflection coefficient, can return a vector or just one"
function reflection_coefficient(ωs::Union{T,AbstractVector{T}}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number

    k_effs = wavenumber(ωs, medium, species; kws...)
    waves = [EffectiveWave(ωs[i], k_effs[i], medium, species; kws...)  for i in eachindex(ωs)]

    return reflection_coefficient(ωs, waves, medium, species; kws...)
end

"The average reflection coefficients"
reflection_coefficient(ωs::AbstractVector{T}, waves::Vector{EffectiveWave{T}}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number =
[ reflection_coefficient(ωs[i], waves[i], medium, species; kws...) for i in eachindex(ωs)]

# reflection_coefficient(ωs::AbstractVector{T}, k_effs::Vector{Complex{T}}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number =
# [ reflection_coefficient(ωs[i], k_effs[i], medium, species; kws...) for i in eachindex(ωs)]

"The average reflection coefficient"
function reflection_coefficient(ω::T, wave_eff::EffectiveWave{T}, medium::Medium{T}, species::Vector{Specie{T}};
        θin::T = 0.0, tol = 1e-8, kws...) where T<:Number

    k = ω/medium.c
    θ_ref = pi - wave_eff.θ_eff - θin
    S = length(species)
    ho = wave_eff.hankel_order

    R = 2.0im/(k*cos(θin)*(k*cos(θin) + wave_eff.k_eff*cos(wave_eff.θ_eff)))
    R = R*sum(
        exp(im*n*θ_ref)*species[l].num_density*wave_eff.amplitudes[n+ho+1,l]*Zn(ω,species[l],medium,n)
    for n=-ho:ho, l=1:S)

    return R
end
