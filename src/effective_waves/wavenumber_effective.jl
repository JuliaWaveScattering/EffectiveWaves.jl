# Here we calculate the effective wavenumber and effective wave amplitudes, without any restriction to the volume fraction of particles and incident wave frequency.

reduce_kvecs(vec::Vector,tol) = vec

function reduce_kvecs(vecs::Vector{Vector{T}},tol::T) where T<:AbstractFloat
    all_inds = collect(eachindex(vecs))
    vecs = map(vecs) do w
        ind_ins = findall([norm(v - w) < T(10)*tol for v in vecs[all_inds]])
        inds = all_inds[ind_ins]
        deleteat!(all_inds,ind_ins)
        isempty(inds) ? [zero(T),-one(T)] :  mean(vecs[inds])
    end
    vecs = deleteat!(vecs, findall([w[2] < -sqrt(tol) for w in vecs]))
end

# include depricated function to find a single effective wavenumber, when in fact there are many. The code is still used in tests and gives many correct results
# include("wavenumber_single.jl")

" Returns all the transmitted effective wavenumbers"
wavenumbers(ω::T, medium::Medium{T}, specie::Specie{T}; kws...) where T<:Number = wavenumbers(ω, medium, [specie]; kws...)

function wavenumbers(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        num_wavenumbers = 8,
        apply_meshing::Bool = true,
        kws...) where T<:Number

    k_effs = wavenumbers_path(ω, medium, species; num_wavenumbers=num_wavenumbers, kws...)
    num_wavenumbers = min(length(k_effs),num_wavenumbers)

    if num_wavenumbers > 1 &&  apply_meshing
        num_wavenumbers = min(length(k_effs),num_wavenumbers)
        k_effs = wavenumbers_mesh(ω, k_effs[1:num_wavenumbers], medium, species; kws...)
    end
    return k_effs
end
