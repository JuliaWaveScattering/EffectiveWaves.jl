"""
    ParticulateMicrostructure

Represents a microstructure filled with multiply species of particles. ParticulateMicrostructure.paircorrelations specifies the pair correlation between each of the species. That is, how the particles are distributed on average.
"""
struct ParticulateMicrostructure{Dim,P<:PhysicalMedium{Dim},Sps<:Species{Dim},PC<:PairCorrelation} <: Microstructure{Dim}
    medium::P
    species::Sps
    paircorrelations::Matrix{PC}
    function ParticulateMicrostructure(medium::P, sps::Sps, ps::AbstractMatrix{PC}) where {Dim, PC <: PairCorrelation, P<:PhysicalMedium{Dim}, Sps <: Species{Dim}}
        #
        if size(ps,1) != length(sps) || size(ps,2) != length(sps)
            @error "the number of rows, and number of columns, of the matrix $paircorrelations needs to be equal to the length of $sps"
        end

        as = [
            s1.separation_ratio * outer_radius(s1) + s2.separation_ratio * outer_radius(s2)
        for s1 in sps, s2 in sps]
        as_pc = [p.minimal_distance for p in ps]

        if !(as ≈ as_pc)
            @warn "The minimal allowed distance between particles defined by the Species is different to that defined by the DiscretePairCorrelation. In this case, the default will be one given by the DiscretePairCorrelation."
        end

        new{Dim,P,Sps,PC}(medium,sps,ps)
    end
end

function Microstructure(medium::PhysicalMedium, sps::Species{Dim}, ps::AbstractMatrix{PC}) where {Dim, PC <: PairCorrelation}
    ParticulateMicrostructure(medium, sps, ps)
end

Microstructure(medium::PhysicalMedium, s::Specie, ps::PairCorrelation; kws...) = Microstructure(medium, [s], [ps][:, :]; kws...)

Microstructure(medium::PhysicalMedium, s::Specie, pc::PairCorrelationType; kws...) = Microstructure(medium, [s], pc; kws...)

function Microstructure(medium::P, sps::Species{Dim}, pc::PairCorrelationType; kws...) where {Dim,P}

    ps = Array{DiscretePairCorrelation}(undef, length(sps), length(sps))

    for i in eachindex(sps), j in eachindex(sps)
        ps[i,j] = if i == j
            DiscretePairCorrelation(sps[i],pc; kws...)
        else
            DiscretePairCorrelation(sps[i],sps[j],pc; kws...)
        end
    end

    return ParticulateMicrostructure(medium,sps,ps)
end

Microstructure(medium::PhysicalMedium, s::Specie; kws...) = Microstructure(medium,[s]; kws...)

"""
    Microstructure(medium::PhysicalMedium, sps::Vector{Specie})

When no pair-correlation is specified for the species, the microstructure will use the default that assumes that particles can not overlap, but, otherwise, their positions are uncorrelated. This is often called \"Hole Correction\"
"""
function Microstructure(medium::PhysicalMedium, sps::Species{Dim}; kws...) where Dim
    ps = [
        DiscretePairCorrelation(s1, s2; kws...)
    for s1 in sps, s2 in sps]

    return ParticulateMicrostructure(medium,sps,ps)
end

import MultipleScattering.get_t_matrices
import MultipleScattering.t_matrix

get_t_matrices(medium::PhysicalMedium, species::Vector{S}, ω::AbstractFloat, Nh::Integer) where S<:Specie = get_t_matrices(medium, [s.particle for s in species], ω, Nh)

t_matrix(s::Specie, medium::PhysicalMedium, ω::AbstractFloat, order::Integer) = t_matrix(s.particle, medium, ω, order)

"""
    Material(region::Shape, micro::Microstructure, [, numberofparticles = Inf])

Creates a material filled with a [`Microstructure`](@ref)'s inside a [`Shape`](@ref).
"""
struct Material{S<:Shape,M<:Microstructure}
    shape::S
    microstructure::M
    numberofparticles::Float64
    # Enforce that the Dims and Types are all the same
    function Material{S,M}(shape::S,micro::M,num::Number = Inf) where {S<:Shape,M<:Microstructure}
        new{S,M}(shape,micro,num)
    end
end

Material(medium::PhysicalMedium,s::Shape,sps::Species) = Material(s, Microstructure(medium,sps))
Material(medium::PhysicalMedium, shape::Shape, specie::Specie) = Material(shape, Microstructure(medium,[specie]))

function Material(shape::S, micro::PM) where {S<:Shape,PM<:ParticulateMicrostructure}

    rmax = maximum(outer_radius.(micro.species))

    # the volume of the shape that contains all the particle centres
    Vn = volume(Shape(shape; addtodimensions = -rmax))

    numberofparticles = sum(
        number_density(s) * Vn
    for s in micro.species)

    Material{S,PM}(shape,micro,numberofparticles)
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
# setupsymmetry(source::AbstractSource, material::Material{S}) where Dim = WithoutSymmetry{Dim}()

# function setupsymmetry(source::PlaneSource{T,3,1}, material::Material{Sphere{T,3}};
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
# function setupsymmetry(psource::PlaneSource{T,Dim}, material::Material{S}) where {T<:AbstractFloat, Dim, S<:Union{Halfspace{T,Dim},Plate{T,Dim}}}
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
