abstract type EffectiveWaveMode{T,Dim} end

"""
    EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T}, wavevector::Array{Complex{T}},amps::Array{Complex{T}})

Is a struct that represents one mode of the effective scattering coefficients for plane wave symmetry.

This wave mode has frequency ω and has has the value ``A e^{i k \\mathbf v \\cdot \\mathbf x }`` at the point ``\\mathbf x``, where ``A`` are the amplitudes, ``\\mathbf v`` is the wavevector and ``k`` is the effective wavenumber.
"""
struct EffectivePlaneWaveMode{T<:AbstractFloat,Dim} <: EffectiveWaveMode{T,Dim}
    ω::T
    basis_order::Int
    wavevector::SVector{Dim,Complex{T}} # the effective wavevector where sum(wavevector.^2) == wavenumber^2
    amplitudes::Array{Complex{T}} # the effective ampliudes
end

function EffectivePlaneWaveMode(ω::T;
        basis_order::Int = 0,
        wavevector::AbstractArray{Complex{T}} = zeros(Complex{T},2),
        amplitudes::AbstractArray{Complex{T}} = zeros(Complex{T},1)
    ) where T<:AbstractFloat

    EffectivePlaneWaveMode(ω, basis_order, SVector(wavevector...), amplitudes)
end

function EffectivePlaneWaveMode(ω::T, wavevector::SVector{Dim,Complex{T}}, amps::Array{Complex{T}}) where {T<:AbstractFloat, Dim}
    EffectivePlaneWaveMode{T,Dim}(ω, Int( (size(amps,1) - 1) / 2 ), wavevector, amps)
end

import Base.zero

zero(W::Type{EffectivePlaneWaveMode{T,Dim}}) where {T<:AbstractFloat,Dim} = EffectivePlaneWaveMode(zero(T),0,zeros(Complex{T},Dim),zeros(Complex{T},1))

"Calculates the effective wavenumbers and return Vector{EffectivePlaneWaveMode}."
# function effective_wavemodes(ω::T, source::AbstractSource{T}, material::Material{Dim,S,Sps}; kws...) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}} # keeps getting "Unreachable reached" error
function effective_wavemodes(ω::T, source::PlaneSource{T,Dim}, material::Material{Dim,S,Sps}; kws...) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}} # without the parametric types we get a "Unreachable reached" error

    k_effs = wavenumbers(ω, source.medium, material.species; kws... )
    wave_effs = [
        effective_wavemode(ω, k_eff, source, material; kws...)
    for k_eff in k_effs]

    return wave_effs
end
