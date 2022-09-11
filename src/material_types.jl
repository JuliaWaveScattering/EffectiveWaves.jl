# """Extract the dimension of the space that this physical property lives in"""
# dim(p::AbstractSymmetry{Dim}) where {Dim} = Dim


"""
    Specie

Represents a set of particles which are all the same. The type of particle is given by `Specie.particle` and the volume fraction this specie occupies is given by `Specie.volume_fraction`.

We can use `Specie.numberofparticles` to specify the number of particles, otherwise for an infinite `Specie.numberofparticles = Inf`.

The minimum distance between any two particles will equal `outer_radius(Specie) * Specie.exclusion_distance`.
"""
struct Specie{Dim,P<:AbstractParticle{Dim}}
    particle::P
    volume_fraction::AbstractFloat
    exclusion_distance::AbstractFloat
end

# Convenience constructor which does not require explicit types/parameters
function Specie(p::AbstractParticle{Dim};
        number_density::AbstractFloat = 0.0,
        volume_fraction::AbstractFloat = number_density * volume(p), exclusion_distance::AbstractFloat = 1.005
    ) where Dim
    Specie{Dim,typeof(p)}(p,volume_fraction,exclusion_distance)
end

function Specie(medium::P,s::S; kws...) where {Dim,P<:PhysicalMedium{Dim},S<:Shape{Dim}}
    Specie(Particle(medium, s); kws...)
end

function Specie(medium::P, radius::AbstractFloat; kws...) where P<:PhysicalMedium
    Specie(Particle(medium, radius); kws...)
end

# Shorthand for all Vectors of species
Species{Dim,P} = Vector{S} where S<:Specie{Dim,P}

"""
    PairCorrelation

Represents the pair correlation between two types of species, which could be the same.
"""
struct PairCorrelation
    "distance between particles centres"
    r::AbstractVector{T} where T <: Number
    "variation of the pair correlation from 1 (uncorrelated case)"
    dp::AbstractVector{T} where T <: Number

    function PairCorrelation(r::AbstractVector,dp::AbstractVector; tol::AbstractFloat = 1e-4)
        if size(dp) != size(r)
            @error "the size of vector of distances `r` (currently $(size(r))) should be the same as the size of the pair-correlation variation `dp` (currently $(size(dp)))."
        end
        if abs(dp[end]) > tol
            @warn "For the pair-correlation to be accurate, we expect it to be long enough (in terms of the distance `r`) such that the particle positions become uncorrelatd. They become uncorrelated when `dp[end]` tends to zero."
        end
        new(r,dp)
    end
end

abstract type AbstractMicrostructure{Dim} end

"""
    ParticulateMicrostructure

Represents potentially multi-species and also holds information on the pair correlation. That is, how the particles are distributed on average.
"""
struct ParticulateMicrostructure{Dim} <: AbstractMicrostructure{Dim}
    species::Species{Dim}
    paircorrelations::AbstractMatrix{PairCorrelation}
    function ParticulateMicrostructure{Dim}(sps::Species{Dim},ps::AbstractMatrix{PairCorrelation}) where Dim
        #
        if size(ps,1) != length(sps) || size(ps,2) != length(sps)
            @error "the number of rows, and number of columns, of the matrix $paircorrelations needs to be equal to the length of $sps"
        end
        new{Dim}(sps,ps)
    end
end

"Returns the volume fraction of the specie."
volume_fraction(s::Specie) = s.volume_fraction
volume_fraction(ss::Species) = sum(volume_fraction.(ss))

import MultipleScattering.volume
volume(s::Specie) = volume(s.particle)

import MultipleScattering.outer_radius

"""
    number_density(s::Specie)

Gives the number of particles per unit volume. Note this is given exactly by `N / V` where `V` is the volume of the region containing the origins of all particles. For consistency, [`volume_fraction`](@ref) is given by `N * volume(s) / V`.
"""
number_density(s::Specie) = s.volume_fraction / volume(s)
# number_density(s::Specie{2}) where {T} = s.volume_fraction / (outer_radius(s.particle)^2 * pi)
# number_density(s::Specie{3}) where {T} = s.volume_fraction / (T(4/3) * outer_radius(s.particle)^3 * pi)
number_density(ss::Species) = sum(number_density.(ss))

outer_radius(s::Specie) = outer_radius(s.particle)

import MultipleScattering.get_t_matrices
import MultipleScattering.t_matrix

get_t_matrices(medium::PhysicalMedium, species::Vector{S}, ω::AbstractFloat, Nh::Integer) where S<:Specie = get_t_matrices(medium, [s.particle for s in species], ω, Nh)

t_matrix(s::Specie, medium::PhysicalMedium, ω::AbstractFloat, order::Integer) = t_matrix(s.particle, medium, ω, order)

"""
    Material(region::Shape, species::Species [, numberofparticles = Inf])

Creates a material filled with [`Specie`](@ref)'s inside `region`.
"""
struct Material{Dim,S<:Shape}
    shape::S
    microstructure::AbstractMicrostructure{Dim}
    numberofparticles::Number
    # Enforce that the Dims and Types are all the same
    function Material{Dim,S}(shape::S,micro::M,num::Number = Inf) where {Dim,S<:Shape{Dim},M<:AbstractMicrostructure{Dim}}
        new{Dim,S}(shape,micro,num)
    end
end

# Convenience constructor which does not require explicit types/parameters
function Material(shape::S,species::Sps) where {Dim,S<:Shape{Dim},Sps<:Species{Dim}}

    @warn "no pair-correlation was specified for the species. Will use the default that particles can not overlap, but are otherwise their positions are uncorrelated. This is often called \"Hole Correction\""

    ps = [
        begin
            a12 = outer_radius(s1) * s1.exclusion_distance + outer_radius(s2) * s2.exclusion_distance
            r = [a12]
            dp = [zero(typeof(a12))]
            PairCorrelation(r,dp)
        end
    for s1 in species, s2 in species]

    micro = ParticulateMicrostructure{Dim}(species,ps)

    return Material(shape,micro)
end

function Material(shape::S,micro::PM) where {Dim,S<:Shape{Dim},PM<:ParticulateMicrostructure{Dim}}

    rmax = maximum(outer_radius.(micro.species))

    # the volume of the shape that contains all the particle centres
    Vn = volume(Shape(shape; addtodimensions = -rmax))

    numberofparticles = round(sum(
        number_density(s) * Vn
    for s in micro.species))

    Material{Dim,S}(shape,micro,numberofparticles)
end

function Material(shape::S,specie::Sp) where {Dim,S<:Shape{Dim},Sp<:Specie{Dim}}
    Material(shape,[specie])
end

import MultipleScattering.PhysicalMedium
import MultipleScattering.Symmetry

PhysicalMedium(s::Specie) = typeof(s.particle.medium)
PhysicalMedium(m::Material) = PhysicalMedium(m.species[1])

Symmetry(m::Material) = Symmetry(m.shape)

# """
#     setupsymmetry(source::AbstractSource, material::Material)
#
# Returns the shared symmetries between the `source` and `materail`.
# """
# setupsymmetry(source::AbstractSource, material::Material{Dim}) where Dim = WithoutSymmetry{Dim}()

# function setupsymmetry(source::PlaneSource{T,3,1}, material::Material{3,Sphere{T,3}};
#         basis_order::Int = 4, ω::T = 0.9) where T
#
#     ls, ms = spherical_harmonics_indices(basis_order)
#
#     gs = regular_spherical_coefficients(source)(basis_order,origin(material.shape),ω)
#     azimuthal_symmetry = norm((ms[n] != 0) ? gs[n] : zero(T) for n in eachindex(gs)) ≈ zero(T)
#
#     return if azimuthal_symmetry
#         AzimuthalSymmetry{3}()
#     else
#         WithoutSymmetry{3}()
#     end
# end
#
# function setupsymmetry(psource::PlaneSource{T,Dim}, material::Material{Dim,S}) where {T<:AbstractFloat, Dim, S<:Union{Halfspace{T,Dim},Plate{T,Dim}}}
#
#     hv = material.shape.normal
#     kv = psource.direction
#
#     if abs(dot(hv, kv)^2) ≈ abs(dot(hv, hv) * dot(kv, kv))
#         # for direct incidence
#         return PlanarAzimuthalSymmetry{Dim}()
#     else
#         return PlanarSymmetry{Dim}()
#     end
# end
