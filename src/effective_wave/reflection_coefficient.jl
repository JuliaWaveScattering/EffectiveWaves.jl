"The average reflection coefficients"
function reflection_coefficients(ωs::Union{T,AbstractVector{T}}, psource::PlaneSource{T,Dim}, material::Material{Dim,Halfspace{T,Dim}}; kws...) where {T<:Number,Dim}

    Rs = map(ωs) do ω
        k_effs = wavenumbers(ω, psource.medium, material.species; kws...)
        w_effs = [
            effective_wavemode(ω, k_eff, psource, material; kws...)
        for k_eff in k_effs]

        reflection_coefficient(ω, w_effs, psource, material; kws...)
    end

    return Rs
end

"Calculates the reflection coefficient from each wave in waves and then sums the results."
reflection_coefficient(ω::T, waves::Vector{EffectivePlaneWaveMode{T}}, psource::PlaneSource{T,Dim}, material::Material{Dim,Halfspace{T,Dim}}; kws...) where {T<:Number,Dim} =
    sum(w -> reflection_coefficient(ω, w, psource, material; kws...), waves)

"Pairs each ω in ωs with each wave in waves to calculate the refleciton coefficients: reflection_coefficient(ω, wave)"
reflection_coefficients(ωs::AbstractVector{T}, waves::Vector{EffectivePlaneWaveMode{T}}, psource::PlaneSource{T,Dim}, material::Material{Dim,Halfspace{T,Dim}}; kws...) where {T<:Number,Dim} =
[ reflection_coefficient(ωs[i], waves[i], psource, material; kws...) for i in eachindex(ωs)]
