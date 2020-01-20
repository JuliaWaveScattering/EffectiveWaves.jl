function wavenumbers_mesh(ω::T, k_effs::Vector{Complex{T}}, medium::PhysicalMedium{T}, species::Species{T};
        dim = 2,
        tol::T = 1e-6,
        mesh_refine::T = T(0.4),
        inner_optimizer = LBFGS(),
        verbose::Bool = false,
        max_Rek::T = maximum(real(k1) for k1 in k_effs),
        min_Rek::T = minimum(real(k1) for k1 in k_effs),
        max_Imk::T = maximum(imag(k1) for k1 in k_effs),
        min_Imk::T = minimum(imag(k1) for k1 in k_effs),
        kws...) where T<:Number

    low_tol = max(1e-4, tol)*minimum( (k1 == k2) ? Inf : abs(k1-k2) for k1 in k_effs, k2 in k_effs) # a tolerance used for a first pass with time_limit

    # the dispersion equation is given by: `dispersion(k1,k2) = 0` where k_eff = k1 + im*k2.
    dispersion = dispersion_equation(ω, medium, species; tol = low_tol, dim=dim, kws...)

    kx = LinRange(min_Rek, max_Rek, Int(round(length(k_effs)/(2*mesh_refine)))) # tree shape makes this division by 2 natural
    ky = LinRange(min_Imk, max_Imk, Int(round(length(k_effs)/(2*mesh_refine))))
    k_mesh = [[x,y] for x in kx, y in ky][:]

    # Include k_effs in the mesh, we do not assume these are solutions
    k_vecs = [[real(keff),imag(keff)] for keff in k_effs]
    k_mesh = [k_mesh; k_vecs]

    # make slightly bigger box to be constrained within
    min_Rek = (min_Rek < zero(T)) ? min_Rek*T(1.1) : min_Rek*T(0.9)
    min_Imk = (min_Imk < zero(T)) ? min_Imk*T(1.1) : min_Imk*T(0.9)
    max_Rek = (max_Rek > zero(T)) ? max_Rek*T(1.001) : max_Rek*T(0.999)
    max_Imk = (max_Imk > zero(T)) ? max_Imk*T(1.001) : max_Imk*T(0.999)
    filter!(keff -> min_Rek <= keff[1] <= max_Rek && min_Imk <= keff[2] <= max_Imk, k_mesh)

    # Find all wavenumbers
    lower = [min_Rek, min_Imk]
    upper = [max_Rek, max_Imk]

    new_ks = [
        optimize(dispersion, lower, upper, kvec,
            Fminbox(inner_optimizer), Optim.Options(x_tol=low_tol, g_tol = low_tol^3)).minimizer
    for kvec in k_mesh]
    deleteat!(new_ks, findall(dispersion.(new_ks) .> low_tol))
    new_ks = reduce_kvecs(new_ks, low_tol)

    filter!(keff -> lower[1] <= keff[1] <= upper[1] && lower[2] <= keff[2] <= upper[2], new_ks)

    # only keep targets which are not already in k_vecs
    # new_ks = [
    #         (findmin([norm(h - kvec) for kvec in k_vecs])[1] > low_tol) ? h : [zero(T),-one(T)]
    # for h in new_ks]
    # new_ks = reduce_kvecs(new_ks, low_tol)

    # Here we refine the new roots
    new_ks = map(new_ks) do k_vec
        # res = optimize(dispersion, k_vec; g_tol = tol^2.0, f_tol = tol^4.0, x_tol=tol)
        res = optimize(dispersion, lower, upper, k_vec, Fminbox(inner_optimizer),
            Optim.Options(g_tol = tol^3.0, x_tol=tol))
        if res.minimum < T(100)*tol || (Optim.converged(res) && res.minimum < T(10)*low_tol)
            res.minimizer
        else
            [zero(T),-one(T)]
        end
    end

    # only keep targets which are not already in k_vecs
    # new_ks = [
    #         (findmin([norm(h - kvec) for kvec in k_vecs])[1] > 10*tol) ? h : [zero(T),-one(T)]
    # for h in new_ks]
    # new_ks = reduce_kvecs(new_ks, T(10)*tol)

    if verbose println("New roots from mesh refiner:",new_ks) end

    # group together wavenumbers which are closer than tol
    k_vecs = reduce_kvecs(new_ks, T(10)*tol)

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

# function detMM!(F,x)
#     F[1] = abs(det(MM(x[1]+im*x[2])))
# end

# Alternative solvers
# res = nlsolve(detMM!,initial_k_eff; iterations = 10000, factor=2.0)
# k_eff_nl = res.zero[1] + im*res.zero[2]
# lower = [0.,0.]; upper = [T(2)*k0,k0]
# result = optimize(dispersion, initial_k0, lower, upper; g_tol = tol^2.0, f_tol = tol^4.0)
