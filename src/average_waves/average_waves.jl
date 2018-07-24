"A type for the ensemble average scattering coefficients.
Here they are discretised in terms of the depth x of the halfspace"
type Scattering_Amplitudes{T<:Real}
    hankel_order::Int # largest hankel order
    x::Vector{T} # spatial mesh
    amplitudes::Array{Complex{T}} # a matrix of the scattering amplitudes, size(A_mat) = (length(x), 2hankel_order +1)
end

Scattering_Amplitudes(M::Int=0, x::AbstractVector{T}=1.0:1.0, as::AbstractArray{Complex{T}}=[1.0+1.0im]) where T<:Number = Scattering_Amplitudes(M,collect(x),collect(as))

function Scattering_Amplitudes(x::AbstractVector{T}, A_mat::Array{Complex{T}}) where T<:Number
    Scattering_Amplitudes(Int((size(A_mat,2)-1)/2), collect(x), A_mat)
end


"returns a function As, where As(k x) is a an array of the ensemble average scattering amplitudes at depth x inside a halfspace.
Note xs are assumed non-dimensional (just as the rest of the package).
As(k x)[m + M + 1,s] is the m-th hankel order and s-th species average scattering coefficient."
function Scattering_Amplitudes(ω::T, xs::AbstractVector{T}, medium::Medium{T}, species::Vector{Specie{T}};
        tol = 1e-8,
        k_eff::Complex{Float64} = wavenumber_low_volfrac(ω, medium, species; tol=tol),
        θin=0.0, kws...) where T<:Number

    As_eff = reduced_amplitudes_effective(ω, k_eff, medium, species;
            θin=θin, tol=tol, kws...)

    k = ω/medium.c
    ho = Int(round(size(As_eff,1)/2-1/2))
    S = length(species)
    θ_eff = transmission_angle(k, k_eff, θin; tol=tol)
    As = [
        im^Float64(m)*exp(-im*m*θ_eff)*As_eff[m+ho+1,s]*exp(im*k_eff*cos(θ_eff)*x)
    for x in xs, m=-ho:ho, s=1:S]

    return Scattering_Amplitudes(ho,xs,As)
end
