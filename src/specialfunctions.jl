"""Define spherical bessel function of the first kind"""
sbesselj(m,x) = sqrt(pi/(2*x)) * besselj(m+1/2,x)

"""Define spherical hankel function of the first kind"""
shankelh1(m,x) = sqrt(pi/(2*x)) * hankelh1(m+1/2,x)

"Derivative of any spherical bessel function"
function diffsbessel(f::Function,n,z)
    return f(n-1,z) - (n+1) * f(n,z) / z
end

"m-th Derivative of any bessel function"
function diffbessel(f::Function,n,z,m::Int)
    if m == 0
        return f(n, z)
    elseif m > 0
        m = m - 1
        return 0.5*(diffbessel(f,n-1,z,m) - diffbessel(f,n+1,z,m))
    else
        error("Can not differentiate a negative number of times")
    end
end

"m-th Derivative of Hankel function of the first kind"
function diffhankelh1(n,z,m::Int)
    if m == 0
        return hankelh1(n, z)
    elseif m > 0
        m = m - 1
        return 0.5*(diffhankelh1(n-1,z,m) - diffhankelh1(n+1,z,m))
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

function kernelN(n::Int,x::Union{T,Complex{T}},y::Union{T,Complex{T}}; dim = 2) where T<:AbstractFloat
    if dim == 2
        h = hankelh1(n,x); dh = diffhankelh1(n,x)
        j = besselj(n,y); dj = diffbesselj(n,y)
    else
        h = shankelh1(n,x); dh = diffsbessel(shankelh1,n,x)
        j = sbesselj(n,y); dj = diffsbessel(sbesselj,n,x)
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
