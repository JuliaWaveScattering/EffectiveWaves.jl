function wavenumbers_mesh(ω::T, k_effs::Vector{Complex{T}}, medium::PhysicalMedium{T,Dim}, species::Species{T};
        symmetry::AbstractSetupSymmetry = PlanarSymmetry{Dim}(),
        tol::T = 1e-5,
        mesh_refine::T = T(0.4),
        inner_optimizer = LBFGS(),
        verbose::Bool = false,
        max_Rek::T = maximum(real(k1) for k1 in k_effs),
        min_Rek::T = minimum(real(k1) for k1 in k_effs),
        max_Imk::T = maximum(imag(k1) for k1 in k_effs),
        min_Imk::T = minimum(imag(k1) for k1 in k_effs),
        kws...) where {T<:Number,Dim}

    eff_medium = effective_medium(medium, species)
    k0 = ω/eff_medium.c
    if isnan(k0) k0 = ω + T(0)*im end
    # abs(k0) can be used to non-dimensionlise k_vec
    kscale = abs(k0)

    low_tol = 1e-2*minimum( (k1 == k2) ? Inf : abs(k1-k2) for k1 in k_effs, k2 in k_effs) / kscale
    low_tol = max(low_tol, 1e-4, tol) # a tolerance used for a first pass with time_limit

    # the dispersion equation is given by: `dispersion(k1,k2) = 0` where k_eff = k1 + im*k2.
    dispersion_dim = dispersion_equation(ω, medium, species, symmetry; tol = low_tol, kws...)
    dispersion(vec::Vector{T}) = dispersion_dim((vec[1] + vec[2] * im) .* kscale)

    # non-dimensionalise
    min_Rek = min_Rek / kscale
    min_Imk = min_Imk / kscale
    max_Rek = max_Rek / kscale
    max_Imk = max_Imk / kscale

    kx = LinRange(min_Rek, max_Rek, Int(round(length(k_effs)/(2*mesh_refine))))# tree shape makes this division by 2 natural
    ky = LinRange(min_Imk, max_Imk, Int(round(length(k_effs)/(2*mesh_refine))))
    k_mesh = [[x,y] for x in kx, y in ky][:]

    # Include k_effs in the mesh, we do not assume these are solutions
    k_vecs = [[real(keff),imag(keff)] for keff in k_effs] ./ kscale
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

    if verbose println("New roots from mesh refiner:",new_ks) end

    # group together wavenumbers which are closer than tol
    k_vecs = reduce_kvecs(new_ks, T(10)*tol)

    # Finally delete unphysical waves, including waves travelling backwards with almost no attenuation. This only is important in the limit of very low frequency or very weak scatterers.
    deleteat!(k_vecs, findall([k_vec[2] < -low_tol for k_vec in k_vecs]))
    deleteat!(k_vecs, findall([-low_tol < abs(k_vec[2])/k_vec[1] < zero(T) for k_vec in k_vecs]))
    # deleteat!(k_vecs, find(k_vec[2] < tol && k_vec[1] < tol for k_vec in k_vecs))

    k_effs = kscale .* [k_vec[1] + k_vec[2]*im for k_vec in k_vecs]
    k_effs = sort(k_effs, by=imag)

    return k_effs
end

# Alternative solvers
# res = nlsolve(detMM!,initial_k_eff; iterations = 10000, factor=2.0)
# k_eff_nl = res.zero[1] + im*res.zero[2]
# lower = [0.,0.]; upper = [T(2)*k0,k0]
# result = optimize(dispersion, initial_k0, lower, upper; g_tol = tol^2.0, f_tol = tol^4.0)
