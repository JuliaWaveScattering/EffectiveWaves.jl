export kernelN2D, kernelN3D # haven't figured out how best to dispath for kernelN
export transmission_angle_wiener, transmission_angle, transmission_direction

"""
    transmission_direction(k_eff::Complex, incident_wavevector::AbstractArray, surface_normal::AbstractArray)

Gives the effective direction `direction` where `sum(direction .^ 2) = 1.0` and the components of `k_eff .* direction` and `incident_wavevector` which are orthogonal to the surface are the same. Note `surface_normal` is the outward pointing normal and the incident medium's wavenumber `k = sqrt(sum(incident_wavevector .^2))`. Below we deduce the result.

Assume that `dot(v,w) = conj(v[i])w[i]` and `surface_normal[i]*surface_normal[i] = 1`. Let `wnp` be orthogonal to `surface_normal`

    wnp =  incident_wavevector -  dot(surface_normal,incident_wavevector) .* surface_normal

Then we know that the effective transmission wavevector has to equal

    k_eff_vec = wnp + α .* surface_normal

where α is determined so that

    sum(x^2 for x in wnp) + α^2  = sum(x^2 for x in k_eff_vec) = k_eff^2
"""
function transmission_direction(k_eff::Complex{T}, incident_wavevector::AbstractArray{CT} where CT <: Union{T,Complex{T}}, surface_normal::AbstractArray{T};
        tol::T = sqrt(eps(T))) where {T<:AbstractFloat,Dim}

    # Let us first make the surface normal a unit vector which is outward pointing
    surface_normal = surface_normal / norm(surface_normal)
    if real(dot(surface_normal,incident_wavevector)) > 0
        surface_normal = -surface_normal
    end

    wnp = incident_wavevector -  dot(surface_normal,incident_wavevector) .* surface_normal
    α = sqrt(k_eff^2 - sum(x^2 for x in wnp))

    # The conditions below guarantee (in order of priority) that either: 1) the wave attenuates when travelling into the material (imag(α) < 0), or 2) the wave does not attenuate, but travels into the material.
    if tol * real(α) > abs(imag(α))
        α = -α # leads to real(α) < 0
    elseif imag(α) > 0
        α = -α # leads to imag(α) < 0
    elseif isnan(α)
        α = -one(T) # covers the cases k_eff = Inf and k_eff = 0.0
        return wnp + α .* surface_normal
    end

    return (wnp + α .* surface_normal) ./ k_eff
end

function transmission_direction(k_eff::Complex{T},  psource::PlaneSource{T,Dim}, material::Material{Dim}; tol::T = sqrt(eps(T))) where {T,Dim}
    transmission_direction(k_eff, psource.direction, material.shape.normal; tol = tol)
end

transmission_angle(pwave::Union{PlaneSource,EffectivePlaneWaveMode}, material::Material) = transmission_angle(pwave, material.shape)

transmission_angle(pwave::Union{PlaneSource,EffectivePlaneWaveMode}, shape::Union{Plate,Halfspace}) = transmission_angle(pwave.direction, shape.normal)

"""
    transmission_angle(wavevector, surface_normal)

Calculates the transmission angle for a wave propagting inside a material which has the outward normal `surface_normal`. The wave propagates in the direction of `wavevector` with a wavenumber given by `k = sqrt(sum(wavevector .^2))`, where `k` should have a positive imaginary part.
"""
transmission_angle(wavevector::AbstractVector,surface_normal::AbstractVector) = transmission_angle(SVector(wavevector...),SVector(surface_normal...))

"""
    transmission_angle(wavevector::StaticVector{3}, surface_normal::StaticVector{3})

Some care is needed when `wavevector` is a complex vector in 3 dimensions.

## The maths

Let's define the normal projection `vn`:

    n = - surface_normal
    vn = dot(n,wavevector) .* n

remembering that in Julia, and mathematics, `dot(v,w) = sum(conj(v[i])*w[i] for i in eachindex(v))`.

We now define the orthogonal component `vo` and its direction `no`:

    vo = wavevector - vn
    no = vo / sqrt(sum(vo.^2))

Note this assumes that `vo` is a real vector, which it should be due to symmetry.

To determine the transmission angle `θT` we calculate the wavenumber

    k = ± sqrt(sum(wavevector .^2))

where we choose the sign in `±` such that `imag(k) >= 0`.

From the above we can write the `wavevector` in the form

    wavevector = k .* (n .* cos(θ) + no .* sin(θ))

which implies that

    k * cos(θ) = dot(conj(n), wavevector)
    k * sin(θ) = dot(conj(no),wavevector) = dot(conj(no),vo) = sqrt(sum(vo.^2))

From these two we can robustly calculate `θT` as

    θT = atan(sqrt(sum(vo.^2)),dot(n,wavevector))

"""
function transmission_angle(wavevector::SVector{3,CT} where CT <: Union{T,Complex{T}}, surface_normal::SVector{3,T}) where {T<:AbstractFloat}
    n = - surface_normal / norm(surface_normal)
    if real(dot(n,wavevector)) < 0
        n = -n
    end

    kcosθ = dot(n,wavevector)

    vo = wavevector - kcosθ .* n
    ksinθ = sqrt(sum(vo.^2))

    return θ = atan(ksinθ,kcosθ)
end


"""
    transmission_angle(wavevector::StaticVector{2}, surface_normal::StaticVector{2})

Specialised to 2 spatial dimensions.
"""
function transmission_angle(wavevector::SVector{2,CT} where CT <: Union{T,Complex{T}}, surface_normal::SVector{2,T}) where {T<:AbstractFloat}
    n = - surface_normal
    no = [-n[2], n[1]] # guarantee θ grows anti-clockwise

    θ = atan(dot(conj(no),wavevector),dot(n,wavevector))

    return θ
end

"""
calculate effective transmission angle θ_eff. We restrict -pi/2 < Re θ_eff < pi/2 for 2 dimensions.
"""
function transmission_angle_wiener(k::Union{T,Complex{T}}, k_eff::Union{T,Complex{T}}, θin) where T<:Number
    # snell(θ::Array{T}) = abs(k*sin(θin) - k_eff*sin(θ[1] + im*θ[2]))
    # result = optimize(snell, [θin,0.]; x_tol= tol, g_tol= tol^2.0)

    θ_eff = asin(k * sin(θin) / k_eff)

    # if !(abs(real(θ_eff - θin)) <= pi/T(2))
    if sign(real(k_eff * cos(θ_eff))) != sign(real(k_eff))
        θ_eff = pi - θ_eff
    end

    return θ_eff
end
