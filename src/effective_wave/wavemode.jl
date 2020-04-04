"""
    transmission_wavevector(k_eff::Complex, incident_wavevector::AbstractArray, surface_normal::AbstractArray)

Gives the effective wavevector `k_eff_vec` where `sum(x^2 for x in k_eff_vec) = k_eff` and the components of `k_eff_vec` and `incident_wavevector` which are orthogonal to the surface are the same. Note `surface_normal` is the outward pointing normal and the incident medium's wavenumber `k = sqrt(sum(incident_wavevector .^2))`. Below we deduce the result.

Assume that `dot(v,w) = conj(v[i])w[i]` and `surface_normal[i]*surface_normal[i] = 1`. Let `wnp` be orthogonal to `surface_normal`

    wnp =  incident_wavevector -  dot(surface_normal,incident_wavevector) .* surface_normal

Then we know that the effective transmission wavevector has to equal

    k_eff_vec = wnp + α .* surface_normal

where α is determined so that

    sum(x^2 for x in wnp) + α^2  = sum(x^2 for x in k_eff_vec) = k_eff^2
"""
function transmission_wavevector(k_eff::Complex{T}, incident_wavevector::AbstractArray{CT} where CT <: Union{T,Complex{T}}, surface_normal::AbstractArray{T};
        tol::T = sqrt(eps(T))) where {T<:AbstractFloat,Dim}

    surface_normal = surface_normal / norm(surface_normal)
    wnp = incident_wavevector -  dot(surface_normal,incident_wavevector) .* surface_normal
    α = sqrt(k_eff^2 - sum(x^2 for x in wnp))

    # The conditions below guarantee (in order of priority) that either: 1) the wave attenuates when travelling into the material (imag(α) < 0), or 2) the wave does not attenuate, but travels into the material.
    if tol * real(α) > abs(imag(α))
        α = -α # leads to real(α) < 0
    elseif imag(α) > 0
        α = -α # leads to imag(α) < 0
    elseif isnan(α)
        α = -one(T) # covers the cases k_eff = Inf and k_eff = 0.0
    end

    return wnp + α .* surface_normal
end

function transmission_wavevector(k_eff::Complex{T},  psource::PlaneSource{T,Dim}, material::Material{Dim}; tol::T = sqrt(eps(T))) where {T,Dim}
    transmission_wavevector(k_eff, psource.direction, material.shape.normal; tol = tol)
end

transmission_angle(pwave::EffectivePlaneWaveMode, shape::Halfspace) = transmission_angle(pwave.wavevector, shape.normal)

transmission_angle(pwave::PlaneSource, shape::Halfspace) = transmission_angle(pwave.direction, shape.normal)

transmission_angle(pwave::Union{PlaneSource,EffectivePlaneWaveMode}, material::Material) = transmission_angle(pwave, material.shape)

transmission_angle(wavevector::Vector,surface_normal::Vector) = transmission_angle(SVector(wavevector...),SVector(surface_normal...))

function transmission_angle(wavevector::SVector{3,CT} where CT <: Union{T,Complex{T}}, surface_normal::SVector{3,T}) where {T<:AbstractFloat}
# Let's define the projection:
    # n = - surface_normal
    # vn = dot(n,wavevector) .* n
# and orthogonal component:
    # vo = wavevector - vn => no = vo / sqrt(sum(vo.^2))
# then need to determine θ_eff such that:
    # k = ± sqrt(sum(wavevector .^2))
    # wavevector = k .* (n .* cos(θ) + no .* sin(θ)) =>
    # k * cos(θ) = dot(conj(n),wavevector)
# and
    # k * sin(θ) = dot(conj(no),wavevector) = dot(conj(no),vo) = sqrt(sum(vo.^2))
    # => θ = atan(sqrt(sum(vo.^2)),dot(n,wavevector))
# where we assume that dot(v,w) = conj(v[i])*w[i]

    n = - surface_normal / norm(surface_normal)
    kcosθ = dot(n,wavevector)

    vo = wavevector - kcosθ .* n
    ksinθ = sqrt(sum(vo.^2))

    return θ = atan(ksinθ,kcosθ)
end

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

"The average effective transmitted scattering amplitudes for a given effective wavenumber k_eff. Assumes there exists only one k_eff.
The function returns an array A, where
AA(x,y,m,s) = im^m*exp(-im*m*θ_eff)*A[m + max_basis_order +1,s]*exp(im*k_eff*(cos(θ_eff)*x + sin(θin)*y))
where (x,y) are coordinates in the halfspace, m-th hankel order, s-th species,  and AA is the ensemble average scattering coefficient."
function mode_amplitudes(ω::T, k_eff::Complex{T}, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, material::Material{Dim,Halfspace{T,Dim}};
        tol::T = 1e-5,
        kws...) where {T<:Number,Dim}

    MM = eigensystem(ω, psource, material; kws...)

    # calculate effective amplitudes
    MM_svd = svd(MM(k_eff))
    if MM_svd.S[end] > tol*T(100)
        @warn("The effective wavenumber gave an eigenvalue $(MM_svd.S[end]), which would be zero if the solution was exact.")
    end

    S = length(material.species)

    A_null = MM_svd.V[:,end] # eignvector of smallest eigenvalue
    A_null = reshape(A_null, (:,S)) # A_null[:,j] = [A[-ho,j],A[-ho+1,j],...,A[ho,j]]

    return A_null
end
