# Here we calculate the effective wavenumber and effective wave amplitudes, without any restriction to the volume fraction of particles and incident wave frequency.

Nn(n::Int,x::Union{T,Complex{T}},y::Union{T,Complex{T}}) where T<:AbstractFloat = x*diffhankelh1(n,x)*besselj(n,y) - y*hankelh1(n,x)*diffbesselj(n,y)

reduce_kvecs(vec::Vector,tol) = vec

function reduce_kvecs(vecs::Vector{Vector{T}},tol::T) where T<:AbstractFloat
    all_inds = collect(eachindex(vecs))
    vecs = map(vecs) do vec
        ind_ins = find(norm(v - vec) < T(10)*tol for v in vecs[all_inds])
        inds = all_inds[ind_ins]
        deleteat!(all_inds,ind_ins)
        isempty(inds) ? [zero(T),-one(T)] :  mean(vecs[inds])
    end
    vecs = deleteat!(vecs, find(vec[2] < -T(10)*tol for vec in vecs))
end

# include depricated function to find a single effective wavenumber, when in fact there are many. The code is still used in tests and gives many correct results
include("wavenumber_single.jl")

" Returns all the transmitted effective wavenumbers"
wavenumbers(ω::T, medium::Medium{T}, specie::Specie{T}; kws...) where T<:Number = wavenumbers(ω, medium, [specie]; kws...)

function wavenumbers(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; apply_meshing::Bool = false, kws...) where T<:Number
    k_effs = wavenumbers_path(ω, medium, species; kws...)
    if length(k_effs) > 1 && apply_meshing
        k_effs = wavenumbers_mesh(ω, k_effs, medium, species; kws...)
    end
    return k_effs
end
