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
function AverageWave(ω::T, medium::Medium{T},specie::Specie{T};
        wave_eff::EffectiveWave{T} = zero(EffectiveWave{T}), kws...) where T<:Number

    k = ω/medium.c
    (X, (MM_quad,b_mat)) = average_wave_system(ω, medium, specie;  kws...);

    M = Int( (size(b_mat,2) - 1)/2 )
    J = length(collect(X)) - 1

    len = (J + 1) * (2M + 1)
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    As = MM_mat\b
    As_mat = reshape(As, (J+1, 2M+1, 1))

    return AverageWave(M, collect(X), As_mat)
end
