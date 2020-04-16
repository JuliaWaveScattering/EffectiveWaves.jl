abstract type AbstractWaveMode{T,Dim} end

abstract type AbstractRegularWaveMode{T,Dim} <: AbstractWaveMode{T,Dim} end

"""
    EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T}, direction::Array{Complex{T}},amps::Array{Complex{T}})

Is a struct that represents one mode of the effective scattering coefficients for plane wave symmetry.

This wave mode has frequency ω and has has the value ``A e^{i k \\mathbf v \\cdot \\mathbf x }`` at the point ``\\mathbf x``, where ``A`` are the amplitudes, ``\\mathbf v`` is the direction and ``k`` is the effective wavenumber.

To recover the average scattering coefficients from a particle `F` use:
    F_{ms} = i^m A[m + max_basis_order +1,s] e^{-i m  θ_eff} * e^{i k \\mathbf v \\cdot \\mathbf x}

where (x,y) are coordinates in the halfspace, m-th hankel order, and s-th species.
"""
struct EffectivePlaneWaveMode{T<:AbstractFloat,Dim} <: AbstractWaveMode{T,Dim}
    ω::T
    wavenumber::Complex{T}
    basis_order::Int
    direction::SVector{Dim,Complex{T}} # the effective direction where sum(direction.^2) == 1
    amplitudes::Array{Complex{T}} # the effective ampliudes
end

function EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T};
        basis_order::Int = 0,
        direction::AbstractArray{Complex{T}} = zeros(Complex{T},2),
        amplitudes::AbstractArray{Complex{T}} = zeros(Complex{T},1)
    ) where T<:AbstractFloat

    EffectivePlaneWaveMode(ω, wavenumber, basis_order, SVector(direction...), amplitudes)
end

function EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T}, direction::SVector{Dim,Complex{T}}, amps::Array{Complex{T}}) where {T<:AbstractFloat, Dim}
    EffectivePlaneWaveMode{T,Dim}(ω, wavenumber, Int( (size(amps,1) - 1) / 2 ), direction, amps)
end

# import Base.zero
#
# zero(W::Type{EffectivePlaneWaveMode{T,Dim}}) where {T<:AbstractFloat,Dim} = EffectivePlaneWaveMode(zero(T), zero(Complex{T}), 0,zeros(Complex{T},Dim),zeros(Complex{T},1))

struct EffectiveAzimuthalWaveMode{T<:AbstractFloat,Dim,P<:PhysicalMedium{T,Dim}} <: AbstractRegularWaveMode{T,Dim}
    ω::T
    wavenumber::Complex{T}
    basis_order::Int
    basis_field_order::Int
    numberofspecies::Int
    amplitudes::Matrix{Complex{T}} # the effective ampliudes, each column is one eigenvector
    function EffectiveAzimuthalWaveMode{T,Dim,P}(ω::T, wavenumber::Complex{T}, basis_order::Int, basis_field_order::Int, numberofspecies::Int, amplitudes::Array{Complex{T}}) where {T,Dim,P<:PhysicalMedium{T,Dim}}

        S = numberofspecies
        L = basis_order
        L1 = basis_order_field

        len = Int(1 - L*(2 + L)*(L - 3*L1 - 2)/3 + L1)

        if size(ampliudes,1) != len
            throw(DimensionMismatch("size(amplitudes,1) does not match the dimensions for a regular eigenvector with azimuthal symmetry."))
        else
            new{T,Dim,P}(ω, wavenumber, basis_order, basis_field_order, numberofspecies, amplitudes)
        end
    end
end


struct EffectiveRegularWaveMode{T<:AbstractFloat,Dim,P<:PhysicalMedium{T,Dim}} <: AbstractRegularWaveMode{T,Dim}
    ω::T
    wavenumber::Complex{T}
    basis_order::Int
    basis_field_order::Int
    numberofspecies::Int
    amplitudes::Matrix{Complex{T}} # the effective ampliudes, each column is one eigenvector
    function EffectiveRegularWaveMode{T,Dim,P}(ω::T, wavenumber::Complex{T}, basis_order::Int, basis_field_order::Int, numberofspecies::Int, amplitudes::Array{Complex{T}}) where {T,Dim,P<:PhysicalMedium{T,Dim}}

        len = (basis_order+1)^2 * (basis_field_order+1)^2 * numberofspecies

        if size(ampliudes,1) != len
            throw(DimensionMismatch("size(amplitudes,1) does not match the dimensions for a regular eigenvector with no symmetry."))
        else
            new{T,Dim,P}(ω, wavenumber, basis_order, basis_field_order, numberofspecies, amplitudes)
        end
    end
end
