"A circular or spherical species of particles with homogenious material properties"
type Specie{T<:Real}
  ρ::T # density
  r::T # radius
  c::Complex{T} # sound speed
  num_density::T # number density
end

"Returns the volume fraction of the specie"
volume_fraction(sp::Specie) = sp.r^2*sp.num_density*pi

function Specie{T}(ρ::T, r::T, c=one(Complex{T}); volfrac::T=0.1*one(T))
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

type Medium{T}
  ρ::T # density
  c::Complex{T} # sound speed
end

Medium(;ρ=1.0, c=1.0+0.0im) = Medium{typeof(ρ)}(ρ,Complex{typeof(ρ)}(c))

Zn{T}(ω::T, p::Specie{T}, med::Medium{T},  m::Int) = Zn(Complex{T}(ω), p, med, m)

"Returns a ratio used in multiple scattering which reflects the material properties of the particles"
function Zn{T}(ω::Complex{T}, p::Specie{T}, med::Medium{T},  m::Int)
    m = T(abs(m))
    ak = p.r*ω/med.c
    # check for material properties that don't make sense or haven't been implemented
    if abs(p.c*p.ρ) == T(NaN)
        error("scattering from a particle with density =$(p.ρ) and phase speed =$(p.c) is not defined")
    elseif abs(med.c*med.ρ) == T(NaN)
        error("wave propagation in a medium with density =$(med.ρ) and phase speed =$(med.c) is not defined")
    elseif abs(med.c) == zero(T)
        error("wave propagation in a medium with phase speed =$(med.c) is not defined")
    elseif abs(med.ρ) == zero(T) && abs(p.c*p.ρ) == zero(T)
        error("scattering in a medium with density $(med.ρ) and a particle with density =$(p.ρ) and phase speed =$(p.c) is not defined")
    end

    # set the scattering strength and type
    if abs(p.c) == T(Inf) || abs(p.ρ) == T(Inf)
        numer = diffbesselj(m, ak)
        denom = diffhankelh1(m, ak)
    elseif abs(med.ρ) == zero(T)
        γ = med.c/p.c #speed ratio
        numer = diffbesselj(m, ak) * besselj(m, γ * ak)
        denom = diffhankelh1(m, ak) * besselj(m, γ * ak)
    else
        q = (p.c*p.ρ)/(med.c*med.ρ) #the impedance
        if q == zero(T)
          numer =  besselj(m, ak)
          denom =  hankelh1(m, ak)
        else
          γ = med.c/p.c #speed ratio
          numer = q * diffbesselj(m, ak)*besselj(m, γ*ak) - besselj(m, ak)*diffbesselj(m, γ*ak)
          denom = q * diffhankelh1(m, ak)*besselj(m, γ*ak) - hankelh1(m, ak)*diffbesselj(m, γ*ak)
        end
    end

    return numer / denom
end

"Derivative of Hankel function of the first kind"
function diffhankelh1(n,z)
  if n!=0
    0.5*(hankelh1(-1 + n, z) - hankelh1(1 + n, z))
  else
    - hankelh1(1, z)
  end
end

"Derivative of Bessel function of first kind"
function diffbesselj(n,z)
  if n!=0
    0.5*(besselj(-1 + n, z) - besselj(1 + n, z))
  else
    - besselj(1, z)
  end
end
