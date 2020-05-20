abstract type AbstractWaveMode{T,Dim} end

abstract type AbstractRegularWaveMode{T,Dim} <: AbstractWaveMode{T,Dim} end


# Each type of symmetry leads to effective eigenvectors with a different length.
eigenvector_length(::PlanarSymmetry{3}; basis_order::Int) =  Int((1 + basis_order)^2)
eigenvector_length(::PlanarSymmetry{2}; basis_order::Int) =  Int((1 + 2*basis_order))

eigenvector_length(::WithoutSymmetry{2}; basis_order::Int, basis_field_order::Int) =  Int( (1 + 2*basis_order) * (1 + 2*basis_field_order))

eigenvector_length(::WithoutSymmetry{3}; basis_order::Int, basis_field_order::Int) =  Int( (1 + basis_order)^2 * (1 + basis_field_order)^2)

eigenvector_length(::AzimuthalSymmetry{3}; basis_order::Int, basis_field_order::Int) =  Int(1 - basis_order*(2 + basis_order)*(basis_order - 3*basis_field_order - 2)/3 + basis_field_order)


"""
    WaveMode(ω::T, wavenumber::Complex{T}, eigenvectors::Array{Complex{T}}, ::SetupSymmetry; kws...)

Returns a concrete subtype of AbstractWaveMode depending on the SetupSymmetry. The returned type should have all the necessary fields to calculate scattered waves (currently not true for EffectivePlanarWaves).
"""
function WaveMode(ω::T, wavenumber::Complex{T}, source::AbstractSource{T}, material::Material{Dim}, eigenvectors::Array{Complex{T}}; kws...) where {T,Dim}

    return EffectiveRegularWaveMode(ω, wavenumber, source, material, eigenvectors; kws...)
end

function WaveMode(ω::T, wavenumber::Complex{T}, psource::PlaneSource{T,Dim,1}, material::Material{Dim,Halfspace{T,Dim}}, eigenvectors::Array{Complex{T}};
    tol::T = 1e-6, kws...) where {T,Dim}

    direction = transmission_direction(wavenumber, (ω / psource.medium.c) * psource.direction, material.shape.normal; tol = tol)

    return EffectivePlaneWaveMode(ω, wavenumber, direction, eigenvectors)
end


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

struct EffectiveRegularWaveMode{T<:AbstractFloat,Dim,P<:PhysicalMedium{T,Dim},S<:AbstractSetupSymmetry{Dim}} <: AbstractRegularWaveMode{T,Dim}
    ω::T
    wavenumber::Complex{T}
    medium::P
    material::Material{Dim}
    eigenvectors::Array{Complex{T}} # the effective eigenvectors, each column is one eigenvector
    basis_order::Int
    basis_field_order::Int
    function EffectiveRegularWaveMode(ω::T, wavenumber::Complex{T}, source::AbstractSource{T}, material::Material{Dim}, eigenvectors::Array{Complex{T}};
        basis_order::Int = 2, basis_field_order::Int = 4
    ) where {T,Dim}

        S = setupsymmetry(source,material)
        P = typeof(source.medium)

        if size(eigenvectors,1) != eigenvector_length(S; basis_order = basis_order, basis_field_order = basis_field_order)
            throw(DimensionMismatch("size(eigenvectors,1) does not match the dimensions for a regular eigenvector with symmetry: $S."))
        elseif size(eigenvectors,2) != length(material.species)
            throw(DimensionMismatch("size(eigenvectors,2) does not match the number of difference species length(material.species) = $(length(material.species)    )."))
        end

        return new{T,Dim,P,typeof(S)}(ω, wavenumber, source.medium, material, eigenvectors, basis_order, basis_field_order)
    end
end
