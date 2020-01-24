"""
    EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T}, wavevector::Array{Complex{T}},amps::Array{Complex{T}})

Is a struct that represents one mode of the effective scattering coefficients for plane wave symmetry.

This wave mode has frequency ω and has has the value ``A e^{i k \\mathbf v \\cdot \\mathbf x }`` at the point ``\\mathbf x``, where ``A`` are the amplitudes, ``\\mathbf v`` is the wavevector and ``k`` is the effective wavenumber.
"""
struct EffectivePlaneWaveMode{T<:AbstractFloat,Dim}
    ω::T
    basis_order::Int
    wavenumber::Complex{T} # the effective wavenumber
    wavevector::SVector{Dim,Complex{T}} # the effective wavevector where sum(wavevector.^2) == 1
    amplitudes::Array{Complex{T}} # the effective ampliudes
end

function EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T};
        basis_order::Int = 0,
        wavevector::Array{Complex{T}} = zeros(Complex{T},2),
        amps::Array{Complex{T}} = zeros(Complex{T},1)) where T<:AbstractFloat

    EffectivePlaneWaveMode(ω, basis_order, wavenumber, wavevector, amps)
end

function EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T}, wavevector::SVector{Dim,Complex{T}},amps::Array{Complex{T}}) where {T<:AbstractFloat, Dim}
    EffectivePlaneWaveMode(ω, Int( (size(amps,1) - 1) / 2 ), wavenumber, wavevector, amps)
end

zero(W::Type{EffectivePlaneWaveMode{T}}) where {T<:AbstractFloat} = EffectivePlaneWaveMode(zero(T),0,[zero(Complex{T})],zero(Complex{T}),zero(Complex{T}))

# effective_wavemodes(ω::T, medium::PhysicalMedium{T}, specie::Specie{T}; kws...) where T<:AbstractFloat =  effective_wavemodes(ω, medium, [specie]; kws...)

"Calculates the effective wavenumbers and return Vector{EffectivePlaneWaveMode}."
function effective_wavemodes(ω::T, source::AbstractSource, material::Material; kws...) where {T<:AbstractFloat,Dim}

    k_effs = wavenumbers(ω, source.medium, material.species; kws... )
    wave_effs = [
        EffectivePlaneWaveMode(ω, k_eff, psource, material; kws...)
    for k_eff in k_effs]

    return wave_effs
end
