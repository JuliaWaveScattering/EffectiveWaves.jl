"A circular or spherical specie of particle with homogenious material properties."
mutable struct Specie{T<:Real}
  ρ::T # density
  r::T # radius
  c::Complex{T} # sound speed
  num_density::T # number density
end

"A homogenious background material."
mutable struct Medium{T}
  ρ::T # density
  c::Complex{T} # sound speed
end

"Returns the volume fraction of the specie."
volume_fraction(sp::Specie) = sp.r^2*sp.num_density*pi

"Returns pressure wave speed, when β is the adiabatic bulk modulus/"
p_speed(ρ,β,μ) = sqrt(β+4*μ/3)/sqrt(ρ)

function Specie(ρ::T, r::T, c=one(Complex{T}); volfrac::T=0.1*one(T)) where T <: AbstractFloat
  Specie{T}(ρ,r,c,volfrac/(T(pi)*r^2.0))
end

function Specie(T=Float64;ρ = 1.0, r = 1.0, c= 1.0+0.0im, volfrac=0.1)
  if T <: Int T = Float64 end
  Specie{T}(T(ρ),T(r),Complex{T}(c),T(volfrac/(pi*r^2.0)))
end
# function Specie(T=Float64; ρ = zero(T), r = one(T), c=one(Complex{T}), volfrac=0.1*one(T))
#   Specie{T}(T(ρ),T(r),c,T(volfrac/(T(pi)*r^2.0)))
# end

# function Specie(T; ρ = zero(Float64), r = one(Float64), c=one(Complex{Float64}), volfrac=0.1*one(Float64))
#   Specie{Float64}(Float64(ρ),Float64(r),c,Float64(volfrac/(Float64(pi)*r^2.0)))
# end

Medium(;ρ=1.0, c=1.0+0.0im) = Medium{typeof(ρ)}(ρ,Complex{typeof(ρ)}(c))
# Medium(;ρ::T = 1.0, c::Union{R,Complex{T}} = 1.0+0.0im) = Medium(ρ,Complex{T}(c))

Nn(n::Int,x::Union{T,Complex{T}},y::Union{T,Complex{T}}) where T<:AbstractFloat =
    x*diffhankelh1(n,x)*besselj(n,y) - y*hankelh1(n,x)*diffbesselj(n,y)
