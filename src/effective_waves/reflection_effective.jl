"The average reflection coefficient, can return a vector or just one"
function reflection_coefficient(ωs::Union{T,AbstractVector{T}}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number

    k_effs = wavenumber(ωs, medium, species; kws...)
    return reflection_coefficient(ωs, k_effs, medium, species; kws...)
end

"The average reflection coefficients"
reflection_coefficient(ωs::AbstractVector{T}, k_effs::Vector{Complex{T}}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number =
[ reflection_coefficient(ωs[i], k_effs[i], medium, species; kws...) for i in eachindex(ωs)]

"The average reflection coefficient"
function reflection_coefficient(ω::T, k_eff::Complex{T}, medium::Medium{T}, species::Vector{Specie{T}};
        θin::T = 0.0, tol = 1e-8, kws...) where T<:Number

    k = ω/medium.c
    θ_eff = transmission_angle(k, k_eff, θin; tol = tol)
    As = reduced_amplitudes_effective(ω, k_eff, medium, species; tol = tol, θin = θin, kws...)
    hankel_order = Int(size(As,1)/2-1/2)

    wave_eff = EffectiveWave(hankel_order, As, k_eff, θ_eff)

    As = As.*scale_amplitudes_effective(ω, wave_eff, medium, species; tol = tol, θin=θin)

    θ_ref = pi - θ_eff - θin
    S = length(species)
    ho = Int((size(As,1)-1)/2) # largest hankel order

    R = 2.0im/(k*cos(θin)*(k*cos(θin) + k_eff*cos(θ_eff)))
    R = R*sum(
        exp(im*n*θ_ref)*species[l].num_density*As[n+ho+1,l]*Zn(ω,species[l],medium,n)
    for n=-ho:ho, l=1:S)

    return R
end
