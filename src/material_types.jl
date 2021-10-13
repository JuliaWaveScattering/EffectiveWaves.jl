export setupsymmetry
export WithoutSymmetry, PlanarSymmetry, PlanarAzimuthalSymmetry, AzimuthalSymmetry, RadialSymmetry

"""
    AbstractSetupSymmetry

An abstract types which dictates the symmetry of the setup. That is, the symmetry shared between the incident wave and the shape of the material.
"""
abstract type AbstractSetupSymmetry{Dim} end

"""
    WithoutSymmetry

A type used to describe materials and incident waves which share no common symmetry. This will lead to the most general regular wave expansion for the eignvectors.
"""
struct WithoutSymmetry{Dim} <: AbstractSetupSymmetry{Dim} end

"""
An incident plane-wave and halfspace material will result in all fields being plane-waves.
"""
abstract type AbstractPlanarSymmetry{Dim} <: AbstractSetupSymmetry{Dim} end

"""
    PlanarSymmetry

A type used to describe materials and incident waves that both share a planar symmetry.
"""
struct PlanarSymmetry{Dim} <: AbstractPlanarSymmetry{Dim} end

"""
For spatial dimension > 2, we can consider problems that have azimuthal symmetry. For example, a plane-wave incident on a sphere.
"""
abstract type AbstractAzimuthalSymmetry{Dim} <: AbstractSetupSymmetry{Dim} end

"""
    AzimuthalSymmetry

A type used to describe materials and incident waves in 3 spatial dimensions that share a symmetry about the azimuth.
"""
struct AzimuthalSymmetry{Dim} <: AbstractAzimuthalSymmetry{Dim} end
AzimuthalSymmetry() = AzimuthalSymmetry{3}()

"""
    RadialSymmetry

A type used to describe materials and incident waves in that are both radially symmetric. That is, the material is a sphere, and the incident wave is radially symmetric around the center of this spherical material.
"""
struct RadialSymmetry{Dim} <: AbstractAzimuthalSymmetry{Dim} end

"""
    PlanarAzimuthalSymmetry

A type used to describe a materail and incident wave which have both [`PlanarSymmetry`](@ref) and [`AzimuthalSymmetry`](@ref).
"""
struct PlanarAzimuthalSymmetry{Dim} <: AbstractPlanarSymmetry{Dim} end
PlanarAzimuthalSymmetry() = PlanarAzimuthalSymmetry{3}()


# """Extract the dimension of the space that this physical property lives in"""
# dim(p::AbstractSetupSymmetry{Dim}) where {Dim} = Dim


"""
    Specie

Represents a set of particles which are all the same. The type of particle is given by `Specie.particle` and the volume fraction this specie occupies is given by `Specie.volume_fraction`.

We can use `Specie.numberofparticles` to specify the number of particles, otherwise for an infinite `Specie.numberofparticles = Inf`.

The minimum distance between any two particles will equal `outer_radius(Specie) * Specie.exclusion_distance`.
"""
struct Specie{T<:AbstractFloat,Dim,P<:AbstractParticle{T,Dim}}
    particle::P
    volume_fraction::T
    exclusion_distance::T
end

# Convenience constructor which does not require explicit types/parameters
function Specie(p::AbstractParticle{T,Dim}; number_density::T = 0.0, volume_fraction::T = number_density * volume(p), exclusion_distance::T = 1.005) where {Dim,T<:AbstractFloat}
    Specie{T,Dim,typeof(p)}(p,volume_fraction,exclusion_distance)
end

function Specie(medium::P,s::S; kws...) where {Dim,T,P<:PhysicalMedium{T,Dim},S<:Shape{T,Dim}}
    Specie(Particle(medium, s); kws...)
end

function Specie(medium::P, radius::T; kws...) where {T,P<:PhysicalMedium{T}}
    Specie(Particle(medium, radius); kws...)
end


# Shorthand for all Vectors of species
Species{T<:AbstractFloat,Dim,P} = Vector{S} where S<:Specie{T,Dim,P}

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
# number_density(s::Specie{T,2}) where {T} = s.volume_fraction / (outer_radius(s.particle)^2 * pi)
# number_density(s::Specie{T,3}) where {T} = s.volume_fraction / (T(4/3) * outer_radius(s.particle)^3 * pi)
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
struct Material{Dim,S<:Shape,Sps<:Species}
    shape::S
    species::Sps
    numberofparticles::Number
    # Enforce that the Dims and Types are all the same
    function Material{Dim,S,Sps}(shape::S,species::Sps,numberofparticles::Number = Inf) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}}
        new{Dim,S,Sps}(shape,species,numberofparticles)
    end
end

# Convenience constructor which does not require explicit types/parameters
function Material(shape::S,species::Sps) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}}

    V = volume(shape)
    numberofparticles = round(sum(
        volume_fraction(s) * V / volume(s)
    for s in species))

    Material{Dim,S,Sps}(shape,species,numberofparticles)
end

function Material(shape::S,specie::Sp) where {T,Dim,S<:Shape{T,Dim},Sp<:Specie{T,Dim}}
    Material(shape,[specie])
end

import MultipleScattering.PhysicalMedium

PhysicalMedium(s::Specie) = typeof(s.particle.medium)
PhysicalMedium(m::Material) = PhysicalMedium(m.species[1])

"""
    setupsymmetry(source::AbstractSource, material::Material)

Returns the shared symmetries between the `source` and `materail`.
"""
setupsymmetry(source::AbstractSource, material::Material{Dim}) where Dim = WithoutSymmetry{Dim}()

function setupsymmetry(source::PlaneSource{T,3,1}, material::Material{3,Sphere{T,3}};
        basis_order::Int = 4, ω::T = 0.9) where T

    ls, ms = spherical_harmonics_indices(basis_order)

    gs = regular_spherical_coefficients(source)(basis_order,origin(material.shape),ω)
    azimuthal_symmetry = norm((ms[n] != 0) ? gs[n] : zero(T) for n in eachindex(gs)) ≈ zero(T)

    return if azimuthal_symmetry
        AzimuthalSymmetry{3}()
    else
        WithoutSymmetry{3}()
    end
end

function setupsymmetry(psource::PlaneSource{T,Dim}, material::Material{Dim,S}) where {T<:AbstractFloat, Dim, S<:Union{Halfspace{T,Dim},Plate{T,Dim}}}

    hv = material.shape.normal
    kv = psource.direction

    if abs(dot(hv, kv)^2) ≈ abs(dot(hv, hv) * dot(kv, kv))
        # for direct incidence
        return PlanarAzimuthalSymmetry{Dim}()
    else
        return PlanarSymmetry{Dim}()
    end
end
