"A type for the effective wave ansatz."
mutable struct EffectiveWave{T<:AbstractFloat}
    hankel_order::Int
    amplitudes::Array{Complex{T}} # the effective ampliudes
    k_eff::Complex{T} # the effective wavenumber
    θ_eff::Complex{T} # the effective transmission angle
end

zero(W::Type{EffectiveWave{T}}) where {T<:AbstractFloat} = EffectiveWave(0,[zero(Complex{T})],zero(Complex{T}),zero(Complex{T}))

effective_waves(ω::T, medium::Medium{T}, specie::Specie{T}; kws...) where T<:AbstractFloat =  effective_waves(ω, medium, [specie]; kws...)

"Calculates the effective wavenumbers and return Vector{EffectiveWave}."
function effective_waves(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; tol::T = 1e-6,
    extinction_rescale::Bool = false, kws...) where T<:AbstractFloat
    # as there will be likely more than 1 k_eff we set extinction to false.

    k_effs = wavenumbers(ω, medium, species; tol = tol, kws... )
    wave_effs = [
        EffectiveWave(ω, k_eff, medium, species; tol = tol, extinction_rescale=extinction_rescale, kws...)
    for k_eff in k_effs]

    return wave_effs
end

function EffectiveWave(ω::T, k_eff::Complex{T}, medium::Medium{T}, species::Vector{Specie{T}};
        θin::T = 0.0, tol::T = 1e-7,
        hankel_order::Int = maximum_hankel_order(ω, medium, species; tol=tol),
        radius_multiplier::T = 1.005,
        method::Symbol = :none,
        extinction_rescale::Bool = true,
        kws...
    ) where T<:AbstractFloat

    k = ω/medium.c
    θ_eff = transmission_angle(k, k_eff, θin; tol = tol)
    if method == :WienerHopf
        if hankel_order > 0
            error("the Wiener Hopf method has only been implemented for monopole scatterers, i.e. hankel order = 0. ")
        end
        amps = wienerhopf_wavevectors(ω, k_eff, medium, species;
            hankel_order=hankel_order, tol = tol,
            radius_multiplier=radius_multiplier,
            θin = θin, kws...
        )
    else
        amps = effective_wavevectors(ω, k_eff, medium, species;
            hankel_order=hankel_order, tol = tol,
            radius_multiplier=radius_multiplier
        )
    end
    if extinction_rescale && method != :WienerHopf
        amps = amps.*scale_amplitudes_effective(ω, wave_eff, medium, species; tol = tol, θin=θin)
    end
    wave_eff = EffectiveWave(hankel_order, amps, k_eff, θ_eff)

    return EffectiveWave(hankel_order, amps, k_eff, θ_eff)
end
