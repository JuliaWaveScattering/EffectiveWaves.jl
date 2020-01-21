"The average reflection coefficients"
function reflection_coefficients(ωs::Union{T,AbstractVector{T}}, medium::PhysicalMedium{T}, species::Species{T}; kws...) where T<:Number

    Rs = map(ωs) do ω
        k_effs = wavenumbers(ω, medium, species; kws...)
        w_effs = [EffectiveWave(ω, k_eff, medium, species; kws...) for k_eff in k_effs]
        reflection_coefficient(ω, w_effs, medium, species; kws...)
    end

    return Rs
end

"Calculates the reflection coefficient from each wave in waves and then sums the results."
reflection_coefficient(ω::T, waves::Vector{EffectiveWave{T}}, medium::PhysicalMedium{T}, species::Species{T}; kws...) where T<:Number =
    sum(w -> reflection_coefficient(ω, w, medium, species; kws...), waves)

"Pairs each ω in ωs with each wave in waves to calculate the refleciton coefficients: reflection_coefficient(ω, wave)"
reflection_coefficients(ωs::AbstractVector{T}, waves::Vector{EffectiveWave{T}}, medium::PhysicalMedium{T}, species::Species{T}; kws...) where T<:Number =
[ reflection_coefficient(ωs[i], waves[i], medium, species; kws...) for i in eachindex(ωs)]
