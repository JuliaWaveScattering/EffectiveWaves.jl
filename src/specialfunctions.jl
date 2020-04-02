export kernelN2D, kernelN3D # haven't figured out how best to dispath for kernelN
export cartesian_to_radial_coordiantes, radial_to_cartesian_coordiantes
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

function kernelN2D(n::Int,x::Union{T,Complex{T}},y::Union{T,Complex{T}}) where T<:AbstractFloat
        h = hankelh1(n,x); dh = diffhankelh1(n,x)
        j = besselj(n,y);  dj = diffbesselj(n,y)

    return x * dh * j - y * h * dj
end

function kernelN3D(n::Int,x::Union{T,Complex{T}},y::Union{T,Complex{T}}) where T<:AbstractFloat
    h = shankelh1(n,x); dh = diffshankelh1(n,x)
    j = sbesselj(n,y);  dj = diffsbesselj(n,y)

    return x * dh * j - y * h * dj
end
