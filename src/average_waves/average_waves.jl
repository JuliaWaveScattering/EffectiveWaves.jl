"A type for the ensemble average scattering coefficients.
Here they are discretised in terms of the depth x of the halfspace"
type AverageWave{T<:Real}
    hankel_order::Int # largest hankel order
    x::Vector{T} # spatial mesh
    amplitudes::Array{Complex{T}} # a matrix of the scattering amplitudes, size(A_mat) = (length(x), 2hankel_order +1)
end

AverageWave(M::Int=0, x::AbstractVector{T}=1.0:1.0, as::AbstractArray{Complex{T}}=[1.0+1.0im]) where T<:Number = AverageWave(M,collect(x),collect(as))

function AverageWave(x::AbstractVector{T}, A_mat::Array{Complex{T}}) where T<:Number
    AverageWave(Int((size(A_mat,2)-1)/2), collect(x), A_mat)
end


"returns a function As, where As(k x) is a an array of the ensemble average scattering amplitudes at depth x inside a halfspace.
Note xs are assumed non-dimensional (just as the rest of the package).
As(k x)[m + M + 1,s] is the m-th hankel order and s-th species average scattering coefficient."
function AverageWave(ω::T, xs::AbstractVector{T}, wave_eff::EffectiveWave{T}, medium::Medium{T}, species::Vector{Specie{T}}) where T<:Number

    k = ω/medium.c
    amps = wave_eff.amplitudes
    ho = wave_eff.hankel_order
    θ_eff = wave_eff.θ_eff

    S = length(species)

    average_amps = [
        im^Float64(m)*exp(-im*m*θ_eff)*amps[m+ho+1,s]*exp(im*wave_eff.k_eff*cos(θ_eff)*x)
    for x in xs, m=-ho:ho, s=1:S]

    return AverageWave(ho,xs,average_amps)
end
