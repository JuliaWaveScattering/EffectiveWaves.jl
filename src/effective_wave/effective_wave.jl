"A type for the effective scattering coefficients for plane wave symmetry."
mutable struct EffectivePlaneWaveMode{T<:AbstractFloat}
    basis_order::Int
    amplitudes::Array{Complex{T}} # the effective ampliudes
    k_eff::Complex{T} # the effective wavenumber
    θ_eff::Complex{T} # the effective transmission angle
end

EffectivePlaneWaveMode(amps::Array{Complex{T}}, k_eff::Complex{T}, θ_eff::Complex{T}) where T<:AbstractFloat = EffectivePlaneWaveMode(Int( (size(amps,1) - 1) / 2 ), amps, k_eff, θ_eff)

zero(W::Type{EffectivePlaneWaveMode{T}}) where {T<:AbstractFloat} = EffectivePlaneWaveMode(0,[zero(Complex{T})],zero(Complex{T}),zero(Complex{T}))

effective_waves(ω::T, medium::PhysicalMedium{T}, specie::Specie{T}; kws...) where T<:AbstractFloat =  effective_waves(ω, medium, [specie]; kws...)

"Calculates the effective wavenumbers and return Vector{EffectivePlaneWaveMode}."
function effective_waves(ω::T, psource::PlaneSource{T,Dim}, material::Material{Dim,Halfspace{T,Dim}}; tol::T = 1e-6, kws...) where {T<:AbstractFloat,Dim}

    k_effs = wavenumbers(ω, psource.medium, material.species; tol = tol, kws... )
    wave_effs = [
        EffectivePlaneWaveMode(ω, k_eff, psource, species; tol = tol, #extinction_rescale=extinction_rescale,
         kws...)
    for k_eff in k_effs]

    return wave_effs
end

function EffectivePlaneWaveMode(ω::T, k_eff::Complex{T}, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, material::Material{Dim,Halfspace{T,Dim}};
        θin::T = 0.0, tol::T = 1e-7,
        method::Symbol = :none,
        extinction_rescale::Bool = true,
        kws...
    ) where {T<:AbstractFloat,Dim}

    k = ω/medium.c
    θ_eff = transmission_angle(k, k_eff, θin; tol = tol)
    if method == :WienerHopf
        amps = wienerhopf_wavemodes(ω, k_eff, medium, species; tol = tol, θin = θin, kws...)
    else
        amps = effective_wavemodes(ω, k_eff, medium, species; tol = tol, kws...)
    end
    wave_eff = EffectivePlaneWaveMode(amps, k_eff, θ_eff)
    if extinction_rescale && method != :WienerHopf
        amps = amps.*scale_amplitudes_effective(ω, wave_eff, medium, species; tol = tol, θin=θin)
    end

    return EffectivePlaneWaveMode(wave_eff.basis_order, amps, k_eff, θ_eff)
end
