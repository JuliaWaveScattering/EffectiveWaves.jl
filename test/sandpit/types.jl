abstract type AbstractParticle{T,Dim} end
abstract type Shape{T<:AbstractFloat,Dim} end

struct Particle{T<:AbstractFloat,Dim,P<:PhysicalProperties,S<:Shape} <: AbstractParticle{T,Dim}
    medium::P
    shape::S
    # Enforce that the Dims and Types are all the same
    function Particle{T,Dim,P,S}(medium::P,shape::S) where {T,Dim,FieldDim,P<:PhysicalProperties{T,Dim,FieldDim},S<:Shape{T,Dim}}
        new{T,Dim,P,S}(medium,shape)
    end
end

struct Specie{T<:Real,Dim,P<:PhysicalProperties,S<:Shape} <: AbstractParticle{T,Dim}
  medium::P
  shape::S
  num_density::T # number density
end


shape(p::Particle) = p.shape
shape(p::CapsuleParticle) = p.outer.shape
