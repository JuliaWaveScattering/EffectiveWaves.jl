"""Define spherical bessel function of the first kind"""
sbesselj(m,x) = sqrt(pi/(2*x)) * besselj(m+1/2,x)

"""Define spherical hankel function of the first kind"""
shankelh1(m,x) = sqrt(pi/(2*x)) * hankelh1(m+1/2,x)

"Derivative of any spherical bessel function"
function diffsbessel(f::Function,n,z)
    return f(n-1,z) - (n+1) * f(n,z) / z
end


# h = 1e-8
# n = 4
# (sbesselj(n, 3.3 + h/2) - sbesselj(n, 3.3 - h/2))/h
# diffsbessel(sbesselj,n,3.3)

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
