type Scattering_Amplitudes{T<:Real}
    hankel_order::Int # largest hankel order
    x::Vector{T} # spatial mesh
    amplitudes::Array{Complex{T}} # a matrix of the scattering amplitudes, size(A_mat) = (length(x), 2hankel_order +1)
end

Scattering_Amplitudes(M::Int=0, x::AbstractVector{T}=1.0:1.0, as::AbstractArray{Complex{T}}=[1.0+1.0im]) where T<:Number = Scattering_Amplitudes(M,collect(x),collect(as))

function Scattering_Amplitudes(x::AbstractVector{T}, A_mat::Array{Complex{T}}) where T<:Number
    Scattering_Amplitudes(Int((size(A_mat,2)-1)/2), collect(x), A_mat)
end
