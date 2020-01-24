export sbesselj, shankelh1, diffsbessel, diffbessel
export gaunt_coefficients, kernelN
export associated_legendre_indices, spherical_harmonics_indices, lm_to_spherical_harmonic_index
export spherical_harmonics, cartesian_to_radial_coordiantes, radial_to_cartesian_coordiantes
export atan


cartesian_to_radial_coordiantes(x::Vector) = cartesian_to_radial_coordiantes(SVector(x...))
radial_to_cartesian_coordiantes(θ::Vector) = radial_to_cartesian_coordiantes(SVector(θ...))

import Base.atan
atan(y::Complex,x::Complex) = - im * log( (x+y*im) / sqrt(x^2+y^2) )

function cartesian_to_radial_coordiantes(x::SVector{3,CT}) where CT
    r = sqrt(sum(x .^2)) # note this should be complex if x is complex
    θ = atan(sqrt(x[1]^2+x[2]^2),x[3])
    φ = atan(x[2], x[1])
    return [r,θ,φ]
end

function cartesian_to_radial_coordiantes(x::SVector{2,CT}) where CT
    r = sqrt(sum(x .^2)) # note this should be complex if x is complex
    θ = atan(x[2], x[1])
    return [r,θ]
end

function radial_to_cartesian_coordiantes(rθφ::SVector{3,CT}) where CT
    r, θ, φ = rθφ
    x = r * sin(θ) * cos(φ)
    y = r * sin(θ) * sin(φ)
    z = r * cos(θ)

    return [x,y,z]
end

function radial_to_cartesian_coordiantes(rθ::SVector{2,CT}) where CT
    r, θ = rθ
    x = r * cos(θ)
    y = r * sin(θ)

    return [x,y]
end

"""
    sbesselj(m,x)

Returns the spherical besselj function. The order is 'm' and the argument is 'x'. Note 'x' can be a complex number.
"""
function sbesselj(m::Int,x::Union{T,Complex{T}}) where T<:AbstractFloat
    if (abs(x) > eps(T))
        return sqrt(pi/(2*x)) * besselj(m+1/2,x)
    else
        return (m > 0 ? zero(typeof(x)) : one(typeof(x)))
    end
end

"""
    shankelh1(m,x)

Returns the spherical hankel function of the first kind. The order is 'm' and the argument is 'x'. Note 'x' can be a complex number.
"""
shankelh1(m,x) = sqrt(pi/(2*x)) * hankelh1(m+1/2,x)

"""
    diffsbessel(f::Function,m,x)

Differentiates the spherical bessel function 'f' (for any spherical bessel). The order is 'm' and the argument is 'x'. Note 'x' can be a complex number.
"""
function diffsbessel(f::Function,n,z)
    return f(n-1,z) - (n+1) * f(n,z) / z
end

"""
    diffbessel(f::Function,m,x,n::Int)

Differentiates 'n' times any bessel function 'f' of order 'm' and at the argument 'x'.
"""
function diffbessel(f::Function,m,z,n::Int)
    if n == 0
        return f(m, z)
    elseif n > 0
        n = n - 1
        return 0.5*(diffbessel(f,m-1,z,n) - diffbessel(f,m+1,z,n))
    else
        error("Can not differentiate a negative number of times")
    end
end

"""
    diffhankelh1(m,x,n::Int)

Differentiates 'n' times the hankelh1 function of order 'm' and at the argument 'x'.
"""
function diffhankelh1(m,z,n::Int)
    if n == 0
        return hankelh1(m, z)
    elseif n > 0
        n = n - 1
        return 0.5*(diffhankelh1(m-1,z,n) - diffhankelh1(m+1,z,n))
    else
        error("Can not differentiate a negative number of times")
    end
end

"Derivative of Hankel function of the first kind"
diffhankelh1(n,z) = 0.5*(hankelh1(-1 + n, z) - hankelh1(1 + n, z))

"m-th Derivative of Hankel function of the first kind"
function diffbesselj(n,z,m::Int)
    if m == 0
        return besselj(n, z)
    elseif m > 0
        m = m - 1
        return 0.5*(diffbesselj(n-1,z,m) - diffbesselj(n+1,z,m))
    else
        error("Can not differentiate a negative number of times")
    end
end

"Derivative of Bessel function of first kind"
diffbesselj(n,z) = 0.5*(besselj(-1 + n, z) - besselj(1 + n, z))

function kernelN(n::Int,x::Union{T,Complex{T}},y::Union{T,Complex{T}}; dim::Int = 2) where T<:AbstractFloat
    if dim == 2
        h = hankelh1(n,x); dh = diffhankelh1(n,x)
        j = besselj(n,y);  dj = diffbesselj(n,y)
    else
        h = shankelh1(n,x); dh = diffsbessel(shankelh1,n,x)
        j = sbesselj(n,y);  dj = diffsbessel(sbesselj,n,y)
    end

    return x * dh * j - y * h * dj
end

"""
    gaunt_coefficients(l1,m1,l2,m2,l3,m3)

A version of the Gaunt coefficients which are used to write the product of two spherical harmonics. If Y_{l,m} is a complex spherical harmonic, with the typical phase conventions from quantum mechanics, then:

    gaunt_coefficients(l1,m1,l2,m2,l3,m3) = 4*π*im^{l2+l3-l1} Integral[Y_{l1,m1}*conj(Y_{l2,m2})*conj(Y_{l3,m3})]

where the integral is over the solid angle.

The most standard gaunt coefficients `G(l1,m1;l2,m2;l3)` are related through the identity:

    4pi * G(l1,m1;l2,m2;l3) = im^(l1-l2-l3) * (-1)^m2 * gaunt_coefficients(l1,m1,l2,-m2,l3,m1+m2)

"""
function gaunt_coefficients(T::Type{<:AbstractFloat},l1::Int,m1::Int,l2::Int,m2::Int,l3::Int,m3::Int)
    # note the wigner3j has only one convention, and is highly symmetric.
    return (one(T)*im)^(l2+l3-l1) * (-T(1))^m1 * sqrt(4pi*(2*l1+1)*(2*l2+1)*(2*l3+1)) *
        wigner3j(T,l1,l2,l3,0,0,0) * wigner3j(T,l1,l2,l3,m1,-m2,-m3)
end
gaunt_coefficients(l1::Int,m1::Int,l2::Int,m2::Int,l3::Int,m3::Int) = gaunt_coefficients(Float64,l1,m1,l2,m2,l3,m3)

lm_to_spherical_harmonic_index(l::Int,m::Int)::Int = l^2 + m + l + 1

function spherical_harmonics_indices(l_max::Int)
    ls = [l for l in 0:l_max for m in -l:l]
    ms = [m for l in 0:l_max for m in -l:l]

    return ls, ms
end

function associated_legendre_indices(l_max::Int)
    ls = [l for l in 0:l_max for m in 0:l]
    ms = [m for l in 0:l_max for m in 0:l]

    return ls, ms
end

"""
`spherical_harmonics(l_max::Int, θ::T, φ::T)`

returns a vector of all spherical harmonics with degree `l <= l_max`. The degree and order (indices) of the elements of the vector are given by `spherical_harmonics_indices(l_max::Int)`.

The associated legendre polynomials are taken from the package GSL.jl.
"""
function spherical_harmonics(l_max::Int, θ::T, φ::T) where T <: AbstractFloat

    ls, ms = associated_legendre_indices(l_max)
    Plm_arr = sf_legendre_array(GSL_SF_LEGENDRE_SPHARM, l_max, cos(θ))[1:length(ls)]

    Ylm_vec = Vector{Complex{T}}(undef, (l_max+1)^2)
    Ylm_vec[1] = Plm_arr[1]

    ind1 = 1
    ind2 = 1
    for i = 1:l_max
        inds1 = (ind1+i):(ind1+2i)
        Ylm_vec[(ind2+i):(ind2+3i)] = [reverse((-1).^ms[inds1[2:end]] .* conj(Plm_arr[inds1[2:end]])); Plm_arr[inds1]]
        ind1 += i
        ind2 += 2i
    end

    ls, ms = spherical_harmonics_indices(l_max)
    Ylm_vec = exp.(ms .* (im*φ)) .* Ylm_vec

    return Ylm_vec
end
