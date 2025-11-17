function wavenumbers_refine(ω::T, medium::Medium, species::Species;
        tol::T = 1e-7,
        k_effs::Vector{Complex{T}} = Complex{T}[],
        kws...) where T<:Number

    # the dispersion equation is given by: `dispersion(k1,k2) = 0` where k_eff = k1 + im*k2.
    dispersion = dispersion_equation(ω, medium, species; tol = tol, kws...)

    # add any specified keffs
    k_vecs = [[real(kp),imag(kp)] for kp in k_effs]

    k_vecs = map(k_vecs) do k_vec
       res = optimize(dispersion, k_vec, Optim.Options(iterations = 10 * Int(abs(round(log(tol)))), g_tol = tol^3.0, x_abstol=tol))
       if res.minimum > T(100)*tol
           [zero(T),-one(T)]
       else
           res.minimizer
       end
    end
    k_vecs = reduce_kvecs(k_vecs, tol)

    # Delete unphysical waves, including waves travelling backwards with almost no attenuation. This only is important in the limit of very low frequency or very weak scatterers.
    deleteat!(k_vecs, findall([k_vec[2] < -sqrt(tol) for k_vec in k_vecs]))

    k_effs = [ k_vec[1] + k_vec[2]*im for k_vec in k_vecs]
    k_effs = sort(k_effs, by=imag)

    return k_effs
end
