function wavenumbers_mesh(ω::T, k_effs::Vector{Complex{T}}, medium::Medium{T}, species::Vector{Specie{T}};
        tol::T = 1e-6,
        hankel_order::Int = maximum_hankel_order(ω, medium, species; tol=tol),
        mesh_refine::T = T(0.4),
        inner_optimizer = LBFGS(),
        verbose::Bool = false,
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; hankel_order = hankel_order),
        max_Rek::T = maximum(real(k1) for k1 in k_effs),
        min_Rek::T = minimum(real(k1) for k1 in k_effs),
        max_Imk::T = maximum(imag(k1) for k1 in k_effs),
        min_Imk::T = minimum(imag(k1) for k1 in k_effs[2:end]),
        kws...) where T<:Number

    k = ω/medium.c
    S = length(species)
    ho = hankel_order

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    function M_component(keff,j,l,m,n)
        (n == m ? 1.0 : 0.0)*(j == l ? 1.0 : 0.0) - 2.0pi*species[l].num_density*t_vecs[l][n+ho+1]*
            kernelN(n-m,k*as[j,l],keff*as[j,l])/(k^2.0-keff^2.0)
    end

    # this matrix is needed to calculate the eigenvectors
    MM(keff::Complex{T}) = reshape(
        [M_component(keff,j,l,m,n) for m in -ho:ho, j = 1:S, n in -ho:ho, l = 1:S]
    , ((2ho+1)*S, (2ho+1)*S))

    detMM2(keff_vec::Array{T}) =  map(x -> real(x*conj(x)), det(MM(keff_vec[1]+im*keff_vec[2])))

    low_tol = max(1e-4, tol)*minimum( (k1 == k2) ? Inf : abs(k1-k2) for k1 in k_effs, k2 in k_effs) # a tolerance used for a first pass with time_limit

    kx = LinRange(min_Rek, max_Rek, Int(round(length(k_effs)/(2*mesh_refine)))) # tree shape makes this division by 2 natural
    ky = LinRange(min_Imk, max_Imk, Int(round(length(k_effs)/(2*mesh_refine))))
    k_mesh = [[x,y] for x in kx, y in ky][:]

    k_vecs = [[real(keff),imag(keff)] for keff in k_effs]

    # make slightly bigger box to be constrained within
    min_Rek = (min_Rek < zero(T)) ? min_Rek*T(1.1) : min_Rek*T(0.9)
    min_Imk = (min_Imk < zero(T)) ? min_Imk*T(1.1) : min_Imk*T(0.9)
    max_Rek = (max_Rek > zero(T)) ? max_Rek*T(1.001) : max_Rek*T(0.999)
    max_Imk = (max_Imk > zero(T)) ? max_Imk*T(1.001) : max_Imk*T(0.999)

    # Find all wavenumbers
    lower = [min_Rek, min_Imk]
    upper = [max_Rek, max_Imk]
    new_ks = [
        optimize(detMM2, lower, upper, kvec,
            Fminbox(inner_optimizer), Optim.Options(x_tol=low_tol, g_tol = low_tol^3)).minimizer
    for kvec in k_mesh]
    deleteat!(new_ks, findall(detMM2.(new_ks) .> low_tol))
    new_ks = reduce_kvecs(new_ks, low_tol/10)

    # only keep targets which are not already in k_vecs
    new_ks = [
            (findmin([norm(h - kvec) for kvec in k_vecs])[1] > low_tol) ? h : [zero(T),-one(T)]
    for h in new_ks]
    new_ks = reduce_kvecs(new_ks, low_tol)

    # Here we refine the new roots
    new_ks = map(new_ks) do k_vec
        # res = optimize(detMM2, k_vec; g_tol = tol^2.0, f_tol = tol^4.0, x_tol=tol)
        res = optimize(detMM2, lower, upper, k_vec, Fminbox(inner_optimizer),
            Optim.Options(g_tol = tol^3.0, x_tol=tol))
        if res.minimum > T(20)*tol
            [zero(T),-one(T)]
        else
            res.minimizer
        end
    end

    new_ks = [
            (findmin([norm(h - kvec) for kvec in k_vecs])[1] > 10*tol) ? h : [zero(T),-one(T)]
    for h in new_ks]
    new_ks = reduce_kvecs(new_ks, tol)

    if verbose println("New roots from mesh refiner:",new_ks) end

    # group together wavenumbers which are closer than tol
    k_vecs = reduce_kvecs([new_ks; k_vecs], tol)

    # Finally delete unphysical waves, including waves travelling backwards with almost no attenuation. This only is important in the limit of very low frequency or very weak scatterers.
    deleteat!(k_vecs, findall([k_vec[2] < -low_tol for k_vec in k_vecs]))
    deleteat!(k_vecs, findall([-low_tol < abs(k_vec[2])/k_vec[1] < zero(T) for k_vec in k_vecs]))
    # deleteat!(k_vecs, find(k_vec[2] < tol && k_vec[1] < tol for k_vec in k_vecs))

    k_effs = [k_vec[1] + k_vec[2]*im for k_vec in k_vecs]
    k_effs = sort(k_effs, by=imag)

    return k_effs
end

# Alternative grouping
# all_inds = collect(eachindex(k_vecs))
# k_vecs = map(k_vecs) do k_vec
#     nk = norm(k_vec)*sqrt(tol)
#     ind_ins = find(norm(ks - k_vec) < nk for ks in k_vecs[all_inds])
#     inds = all_inds[ind_ins]
#     deleteat!(all_inds,ind_ins)
#     isempty(inds) ? [zero(T),-one(T)] :  mean(k_vecs[inds])
# end

# alternative grouping methods

    # digs = 3 - Int(round(log(10,sqrt(tol)))) # number of digits to round to check if equal
    # k_vecs = map(groupby(k_vec -> round.(k_vec,digs), k_vecs)) do g
    #     mean(g)
    # end

# detMM2(keff_vec::Array{T}) =  map(x -> real(x*conj(x)), det(MM(keff_vec[1]+im*keff_vec[2])))

# function detMM!(F,x)
#     F[1] = abs(det(MM(x[1]+im*x[2])))
# end

# Alternative solvers
# res = nlsolve(detMM!,initial_k_eff; iterations = 10000, factor=2.0)
# k_eff_nl = res.zero[1] + im*res.zero[2]
# lower = [0.,0.]; upper = [T(2)*k0,k0]
# result = optimize(detMM2, initial_k0, lower, upper; g_tol = tol^2.0, f_tol = tol^4.0)
