"A type for the effective wave ansatz."
mutable struct EffectiveWave{T<:AbstractFloat}
    basis_order::Int
    amplitudes::Array{Complex{T}} # the effective ampliudes
    k_eff::Complex{T} # the effective wavenumber
    θ_eff::Complex{T} # the effective transmission angle
end

EffectiveWave(amps::Array{Complex{T}}, k_eff::Complex{T}, θ_eff::Complex{T}) where T<:AbstractFloat = EffectiveWave(Int( (size(wave_eff.amplitudes,1) - 1) / 2 ), amps, k_eff, θ_eff)

zero(W::Type{EffectiveWave{T}}) where {T<:AbstractFloat} = EffectiveWave(0,[zero(Complex{T})],zero(Complex{T}),zero(Complex{T}))

effective_waves(ω::T, medium::PhysicalMedium{T}, specie::Specie{T}; kws...) where T<:AbstractFloat =  effective_waves(ω, medium, [specie]; kws...)

"Calculates the effective wavenumbers and return Vector{EffectiveWave}."
function effective_waves(ω::T, medium::Acoustic{T,2}, species::Vector{Specie{T}}; tol::T = 1e-6,
    extinction_rescale::Bool = false, kws...) where T<:AbstractFloat
    # as there will be likely more than 1 k_eff we set extinction to false.

    k_effs = wavenumbers(ω, medium, species; tol = tol, kws... )
    wave_effs = [
        EffectiveWave(ω, k_eff, medium, species; tol = tol, extinction_rescale=extinction_rescale, kws...)
    for k_eff in k_effs]

    return wave_effs
end

function EffectiveWave(ω::T, k_eff::Complex{T}, medium::Acoustic{T,2}, species::Vector{Specie{T,2}};
        θin::T = 0.0, tol::T = 1e-7,
        # basis_order::Int = 2, #maximum_basis_order(ω, medium, species; tol=tol),
        # radius_multiplier::T = 1.005,
        method::Symbol = :none,
        extinction_rescale::Bool = true,
        kws...
    ) where T<:AbstractFloat

    k = ω/medium.c
    θ_eff = transmission_angle(k, k_eff, θin; tol = tol)
    if method == :WienerHopf
        amps = wienerhopf_wavevectors(ω, k_eff, medium, species; tol = tol, θin = θin, kws...)
    else
        amps = effective_wavevectors(ω, k_eff, medium, species; tol = tol, kws...)
    end
    wave_eff = EffectiveWave(amps, k_eff, θ_eff)
    if extinction_rescale && method != :WienerHopf
        amps = amps.*scale_amplitudes_effective(ω, wave_eff, medium, species; tol = tol, θin=θin)
    end

    return EffectiveWave(wave_eff.basis_order, amps, k_eff, θ_eff)
end
