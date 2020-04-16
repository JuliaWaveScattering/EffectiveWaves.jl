abstract type EffectiveWaveMode{T,Dim} end

"""
    EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T}, wavevector::Array{Complex{T}},amps::Array{Complex{T}})

Is a struct that represents one mode of the effective scattering coefficients for plane wave symmetry.

This wave mode has frequency ω and has has the value ``A e^{i k \\mathbf v \\cdot \\mathbf x }`` at the point ``\\mathbf x``, where ``A`` are the amplitudes, ``\\mathbf v`` is the direction and ``k`` is the effective wavenumber.
"""
struct EffectivePlaneWaveMode{T<:AbstractFloat,Dim} <: EffectiveWaveMode{T,Dim}
    ω::T
    wavenumber::Complex{T}
    basis_order::Int
    direction::SVector{Dim,Complex{T}} # the effective direction where sum(direction.^2) == 1
    amplitudes::Array{Complex{T}} # the effective ampliudes
end

function EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T};
        basis_order::Int = 0,
        wavevector::AbstractArray{Complex{T}} = zeros(Complex{T},2),
        amplitudes::AbstractArray{Complex{T}} = zeros(Complex{T},1)
    ) where T<:AbstractFloat

    EffectivePlaneWaveMode(ω, wavenumber, basis_order, SVector(wavevector...), amplitudes)
end

function EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T}, wavevector::SVector{Dim,Complex{T}}, amps::Array{Complex{T}}) where {T<:AbstractFloat, Dim}
    EffectivePlaneWaveMode{T,Dim}(ω, wavenumber, Int( (size(amps,1) - 1) / 2 ), wavevector, amps)
end

import Base.zero

zero(W::Type{EffectivePlaneWaveMode{T,Dim}}) where {T<:AbstractFloat,Dim} = EffectivePlaneWaveMode(zero(T), zero(Complex{T}), 0,zeros(Complex{T},Dim),zeros(Complex{T},1))
