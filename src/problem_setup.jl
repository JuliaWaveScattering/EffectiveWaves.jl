# I need to better understand the general case before designing a nice dispatch system...

"""
An abstract types which dictates the problem to be solved. For example, all problem with two spatial dimensions, all of which use 2D effective wavenumbers, we have the type TwoDimensions{T} <: PhysicalSetup{T}. In particular, for planar interface, we then have PlaneWaves2D{T} <: TwoDimensions{T}.
"""
abstract type ProblemSetup{T<:AbstractFloat, Dim<:Integer} end

abstract type PlaneWaves3D{T} <: ProblemSetup{T,3} end


"""Extract the dimension of the space that this physical property lives in"""
dim(p::ProblemSetup{T,Dim}) where {Dim,T} = Dim

"""
A setup with an incident plane wave, from the left, on a plate region filled with particles.
"""
struct Plate2D{T} <: ProblemSetup{T,2}
    """The plate is the region defined by x in the interval [x1, x2]"""
    x1::T
    x2::T
    """The background medium"""
    medium::PhysicalMedium{T}
    """the species of the particles"""
    species::Vector{Specie{T}}
    """maximum order of hankel functions"""
    basis_order::Int
    """minimal particle distance multiplier"""
    radius_multiplier::T
    """orientation of wave-vector (k*cos(θin), k*sin(θin))"""
    θin::T
end

"""
A setup with an incident plane wave, from the left, on a plate region filled with particles.
"""
struct Sphere{T} <: ProblemSetup{T,3}
    centre::T
    radius::T
    """The background medium"""
    medium::PhysicalMedium{T}
    """the species of the particles"""
    species::Vector{Specie{T}}
    """maximum order of hankel functions"""
    basis_order::Int
    """minimal particle distance multiplier"""
    radius_multiplier::T
    """orientation of wave-vector (k*cos(θin), k*sin(θin))"""
    θin::T
end
