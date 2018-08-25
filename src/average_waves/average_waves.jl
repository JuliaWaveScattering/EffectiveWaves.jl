"A type for the ensemble average scattering coefficients.
Here they are discretised in terms of the depth x of the halfspace"
type AverageWave{T<:AbstractFloat}
    hankel_order::Int # largest hankel order
    X::Vector{T} # spatial mesh
    amplitudes::Array{Complex{T}} # a matrix of the scattering amplitudes, size(A_mat) = (length(x), 2hankel_order +1)
    # Enforce that the dimensions are correct
    function AverageWave{T}(hankel_order::Int, X::Vector{T}, amplitudes::Array{Complex{T}}) where T <: AbstractFloat
        if (length(X), 2*hankel_order+1) != size(amplitudes)[1:2]
            error("The amplitudes of AverageWave does not satisfy size(amplitudes)[1:2] == (length(X), 2*hankel_order+1)")
        end
        new(hankel_order,X,amplitudes)
    end
end

AverageWave(M::Int, X::AbstractVector{T}, as::AbstractArray{Complex{T}}) where T<:AbstractFloat = AverageWave{T}(M,collect(X),collect(as))

function AverageWave(X::AbstractVector{T}, A_mat::Array{Complex{T}}) where T<:Number
    AverageWave(Int((size(A_mat,2)-1)/2), collect(x), A_mat)
end


"Calculates an AverageWave from one EffectiveWave"
function AverageWave(k::T, wave_eff::EffectiveWave{T}, X::AbstractVector{T}, X_match::T = zero(T)) where T<:Number

    amps = wave_eff.amplitudes
    ho = wave_eff.hankel_order
    θ_eff = wave_eff.θ_eff

    S = size(amps,2)

    average_amps = [
        im^T(m)*exp(-im*m*θ_eff)*amps[m+ho+1,s]*exp(im*wave_eff.k_eff*cos(θ_eff)*x/k)
    for x in X, m=-ho:ho, s=1:S]

    return AverageWave(ho,X,average_amps)
end

"Numerically solved the integral equation governing the average wave. Optionally can use wave_eff to approximate the wave away from the boundary."
function AverageWave(ω::T, medium::Medium{T}, specie::Specie{T};
        radius_multiplier::T = 1.005,
        X::T = [zero(T)], tol::T = T(1e-5),
        wave_effs::Vector{EffectiveWave{T}} = [zero(EffectiveWave{T})], kws...) where T<:Number

    k = real(ω/medium.c)
    a12k = T(2)*radius_multiplier*specie.r

    if X == [zero(T)]
        if maximum(abs(w.k_eff) for w in wave_effs) == zero(T)
            error("Either provide the mesh X = .. or the effective waves wave_effs = .., from which we can estimate a mesh X")
        else L, X =  X_match_waves(k, wave_effs, a12k; tol = tol)
        end
    end
    (MM_quad,b_mat) = average_wave_system(ω, X, medium, specie;  kws...);

    M = Int( (size(b_mat,2) - 1)/2 )
    J = length(collect(X)) - 1

    len = (J + 1) * (2M + 1)
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    As = MM_mat\b
    As_mat = reshape(As, (J+1, 2M+1, 1))

    return AverageWave(M, collect(X), As_mat)
end

"note that this uses the non-dimensional X = k*depth"
function average_wave_system(ω::T, X::AbstractVector{T}, medium::Medium{T}, specie::Specie{T};
        θin::Float64 = 0.0, tol::T = 1e-6,
        radius_multiplier::T = 1.005,
        hankel_order::Int = maximum_hankel_order(ω, medium, [specie]; tol = 1000*tol)
    ) where T<:AbstractFloat

    k = real(ω/medium.c)
    a12k = radius_multiplier*T(2)*real(k*specie.r);
    M = hankel_order;

    J = length(collect(X))
    h = X[2] - X[1]

    Z = OffsetArray{Complex{Float64}}(-M:M);
    for m = 0:M
        Z[m] = Zn(ω,specie,medium,m)
        Z[-m] = Z[m]
    end

    σ =  trap_scheme(X) # integration scheme: trapezoidal
    PQ_quad = intergrand_kernel(X, a12k; θin = θin, M = M);

    MM_quad = [
        specie.num_density*Z[n]*σ[j]*PQ_quad[l,m+M+1,j,n+M+1] + k^2*( (m==n && j==l) ? 1.0+0.0im : 0.0+0.0im)
    for  l=1:J, m=-M:M, j=1:J, n=-M:M];

    b_mat = [ -k^2*exp(im*X[l]*cos(θin))*exp(im*m*(pi/2.0 - θin)) for l = 1:J, m = -M:M]

    return (MM_quad,b_mat)
end
