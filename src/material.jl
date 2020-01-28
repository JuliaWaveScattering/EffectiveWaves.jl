
"Represents a set of particles."
struct Specie{T<:AbstractFloat,Dim,P<:AbstractParticle{T,Dim}}
    particle::P
    volume_fraction::T
    numberofparticles::Int
    exclusion_distance::T
end

# Convenience constructor which does not require explicit types/parameters
function Specie(p::AbstractParticle{T,Dim}; volume_fraction::T = 0.0, exclusion_distance::T = 1.005, numberofparticles::Int = -1) where {Dim,T<:AbstractFloat}
    Specie{T,Dim,typeof(p)}(p,volume_fraction,numberofparticles,exclusion_distance)
end

# Shorthand for all Vectors of species
Species{T<:AbstractFloat,Dim,P} = Vector{S} where S<:Specie{T,Dim,P}

"Returns the volume fraction of the specie."
volume_fraction(s::Specie) = s.volume_fraction
volume_fraction(ss::Species) = sum(volume_fraction.(ss))

import MultipleScattering.outer_radius

"Returns the number density of the specie."
number_density(s::Specie{T,2}) where {T} = s.volume_fraction / (outer_radius(s.particle)^2 * pi)
number_density(s::Specie{T,3}) where {T} = s.volume_fraction / (T(4/3) * outer_radius(s.particle)^3 * pi)
number_density(ss::Species) = sum(number_density.(ss))


outer_radius(s::Specie) = outer_radius(s.particle)

import MultipleScattering.get_t_matrices
import MultipleScattering.t_matrix

get_t_matrices(medium::PhysicalMedium, species::Vector{S}, ω::AbstractFloat, Nh::Integer) where S<:Specie = get_t_matrices(medium, [s.particle for s in species], ω, Nh)

t_matrix(s::Specie, medium::PhysicalMedium, ω::AbstractFloat, order::Integer) = t_matrix(s.particle, medium, ω, order)

"""
    Material(region::Shape, species::Species)

Creates a material filled with [`Specie`](@ref)'s inside a region.
"""
struct Material{Dim,S<:Shape,Sps<:Species}
    shape::S
    species::Sps
    # Enforce that the Dims and Types are all the same
    function Material{Dim,S,Sps}(shape::S,species::Sps) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}}
        new{Dim,S,Sps}(shape,species)
    end
end

# Convenience constructor which does not require explicit types/parameters
function Material(shape::S,species::Sps) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}}
    Material{Dim,S,Sps}(shape,species)
end

function Material(shape::S,specie::Sp) where {T,Dim,S<:Shape{T,Dim},Sp<:Specie{T,Dim}}
    Material{Dim,S,Vector{Sp}}(shape,[specie])
end
