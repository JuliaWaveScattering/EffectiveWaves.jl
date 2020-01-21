reduce_kvecs(vec::Vector,tol) = vec

function reduce_kvecs(vecs::Vector{Vector{T}},tol::T) where T<:AbstractFloat
    all_inds = collect(eachindex(vecs))
    vecs = map(vecs) do w
        ind_ins = findall([norm(v - w) < tol for v in vecs[all_inds]])
        inds = all_inds[ind_ins]
        deleteat!(all_inds,ind_ins)
        isempty(inds) ? [zero(T),-one(T)] :  mean(vecs[inds])
    end
    vecs = deleteat!(vecs, findall([w[2] < -sqrt(tol) for w in vecs]))

    return vecs
end
