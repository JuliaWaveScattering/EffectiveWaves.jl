abstract type AbstractWaveMode{T,Dim} end

abstract type AbstractRegularWaveMode{T,Dim} <: AbstractWaveMode{T,Dim} end


# Each type of symmetry leads to effective eigenvectors with a different length.
eigenvector_length(::PlanarSymmetry{3}; basis_order::Int) =  Int((1 + basis_order)^2)
eigenvector_length(::PlanarSymmetry{2}; basis_order::Int) =  Int((1 + 2*basis_order))

eigenvector_length(::WithoutSymmetry{2}; basis_order::Int, basis_field_order::Int) =  Int( (1 + 2*basis_order) * (1 + 2*basis_field_order))

eigenvector_length(::WithoutSymmetry{3}; basis_order::Int, basis_field_order::Int) =  Int( (1 + basis_order)^2 * (1 + basis_field_order)^2)

eigenvector_length(::AzimuthalSymmetry{3}; basis_order::Int, basis_field_order::Int) =  Int(1 - basis_order*(2 + basis_order)*(basis_order - 3*basis_field_order - 2)/3 + basis_field_order)


"""
    EffectivePlaneWaveMode{T<:AbstractFloat,Dim} <: AbstractWaveMode{T,Dim}

Is a struct that represents a mode of the average scattering coefficients for plane wave symmetry. That is, when using an incident plane wave and a plate or halfspace for the material geometry.

This wave mode has frequency ``\\omega`` and has has the value ``A e^{i k \\mathbf v \\cdot \\mathbf x }`` at the point ``\\mathbf x``, where ``A`` are the eigenvectors, ``\\mathbf v`` is the direction and ``k`` is the effective wavenumber. We represent these quantities as

  * ``\\omega = `` `EffectivePlaneWaveMode.ω`
  * ``k = `` `EffectivePlaneWaveMode.wavenumber`
  * ``v = `` `EffectivePlaneWaveMode.direction`
  * ``A = `` `EffectivePlaneWaveMode.eigenvectors`
  * ``M = `` `EffectivePlaneWaveMode.basis_order`

To recover the average scattering coefficients ``F`` from a particle use:

``F_{ns} = i^m A_{n s} e^{-i n  θ} ``

where ``n`` is multi-index for 3 dimensions and ``s`` ranges from 1 to the number of species.

In 2D, inconvieniently, an extra factor of ``e^{i k \\mathbf v \\cdot \\mathbf x}`` needs to be multiplied on the right of the above equation where `θ` is calculated from [`transmission_angle`](@ref).
"""
struct EffectivePlaneWaveMode{T<:AbstractFloat,Dim} <: AbstractWaveMode{T,Dim} #add: P<:PhysicalMedium{T,Dim}
    ω::T
    wavenumber::Complex{T}
    basis_order::Int
    direction::SVector{Dim,Complex{T}} # the effective direction where sum(direction.^2) == 1
    eigenvectors::Array{Complex{T}} # the effective ampliudes
end

function EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T};
        basis_order::Int = 0,
        direction::AbstractArray{Complex{T}} = zeros(Complex{T},2),
        eigenvectors::AbstractArray{Complex{T}} = zeros(Complex{T},1)
    ) where T<:AbstractFloat

    EffectivePlaneWaveMode(ω, wavenumber, basis_order, SVector(direction...), eigenvectors)
end

"""
    EffectivePlaneWaveMode(ω, wavenumber::Complex, direction::Array{Complex}, eigenvectors::Array{Complex})

A convienient way to form the type `EffectivePlaneWaveMode`.
"""
function EffectivePlaneWaveMode(ω::T, wavenumber::Complex{T}, direction::SVector{Dim,Complex{T}}, eigvectors::Array{Complex{T}}) where {T<:AbstractFloat, Dim}
    EffectivePlaneWaveMode{T,Dim}(ω, wavenumber, basislength_to_basisorder(Acoustic{T,Dim},size(eigvectors,1)), direction, eigvectors)
end


"""
    EffectiveRegularWaveMode{T,Dim,P<:PhysicalMedium,S<:AbstractSetupSymmetry} <: AbstractRegularWaveMode{T,Dim}

Is a struct that represents a mode of the average scattering coefficients ``F`` in the form

``F_{n s} (\\mathbf r) =  \\sum_{n_1 p} A_{p n n_1 s} \\mathrm v_{n_1}(k \\mathbf r)``

where ``\\mathbf r`` is a position vector in space, ``k`` the wavenumber, ``p`` is used to count the number of eigenvectors, ``s`` the number of species, ``n`` and ``n_1`` are multi-indices, and ``\\mathrm v_{n_1}`` is a regular spherical wave of order ``n_1`` which in 3D is given by

``\\mathrm v_{n_1}(\\mathbf r) = \\mathrm j_{\\ell_1} (k r) \\mathrm Y_{\\ell_1 m_1}(\\hat{\\mathbf r})``

and ``\\mathrm Y_{\\ell_1 m_1}`` is a spherical harmonic, ``n_1 = (\\ell_1,m_1)`` where we always have that ``\\ell_1 \\geq 0`` and ``\\ell_1 \\geq |m_1|`` according to the conventions of [spherical harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics).

To translate the mathematics to Julia code we use

  * ``\\omega = `` `EffectiveRegularWaveMode.ω`
  * ``k = `` `EffectiveRegularWaveMode.wavenumber`
  * ``A = `` `EffectiveRegularWaveMode.eigenvectors`
  * ``n = `` `EffectiveRegularWaveMode.basis_order`
  * ``n_1 = `` `EffectiveRegularWaveMode.basis_field_order`

"""
struct EffectiveRegularWaveMode{T<:AbstractFloat,Dim,P<:PhysicalMedium{T,Dim},S<:AbstractSetupSymmetry{Dim}} <: AbstractRegularWaveMode{T,Dim}
    ω::T
    wavenumber::Complex{T}
    medium::P
    material::Material{Dim}
    eigenvectors::Array{Complex{T}} # the effective eigenvectors, each column is one eigenvector
    basis_order::Int
    basis_field_order::Int
    function EffectiveRegularWaveMode(ω::T, wavenumber::Complex{T}, source::AbstractSource{T}, material::Material{Dim}, eigenvectors::Array{Complex{T}};
        basis_order::Int = 2, basis_field_order::Int = 4, kws...
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
