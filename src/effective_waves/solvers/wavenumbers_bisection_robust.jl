# NOTE: PlanarAzimuthalSymmetry() does not included all possible wavenumbers
function wavenumbers_bisection_robust(ωs::T, micro::Microstructure{Dim};
        symmetry::AbstractSymmetry{Dim} = PlanarAzimuthalSymmetry{Dim}(),
        tol::T = 1e-5,
        num_wavenumbers = 1,
        box_k::Vector{Vector{T}} = box_keff(ω, micro; num_wavenumbers = num_wavenumbers),
        # bisection_mesh_points::Int = 2 * Int(round(- log(tol))) + num_wavenumbers,
        bisection_iteration::Int = 2,
        fixedparameters::Optim.FixedParameters = NelderMeadparameters(),
        optimoptions::Optim.Options{T} = Optim.Options(
            g_tol = tol^T(3), x_abstol=tol^T(2)),
        kws...) where {T,Dim}

    # exact low-frequency effective wavenumber is known
    # ko = real(ω / effective_medium(micro).c)

    kφ = wavenumber_low_volumefraction(ω, micro;
        verbose = false
    )

    disp = dispersion_complex(ω, micro, symmetry; basis_order = basis_order)

    freal(x, y)::T = real(disp(x + y*im))
    fimag(x, y)::T = imag(disp(x + y*im))

    x = box_k[1][1]:real(kφ):box_k[1][2];
    y = (box_k[2][1] - imag(kφ)):imag(kφ):box_k[2][2];

    # axes=[range ω, real part of kz, imag part of kz]
    # axes = [1:3,0:4,-4:4]

    # Doesn't work as real part has tiny loop around the principal root.
        # #Intersection of both real and imag surface of dispersion equation
        #     axes = [Axis(ω_mesh,"ω"), Axis(x,"x"),Axis(y,"y")]
        #     f_combined(x...) = (freal(x...), fimag(x...))
        #     Intersectmdbm = MDBM_Problem(f_combined,axes)

    axes = [Axis(x,"x"),Axis(y,"y")]
    Intersectmdbm = MDBM_Problem(fimag,axes)

    solve!(Intersectmdbm,bisection_iteration)
    x_sol, y_sol = getinterpolatedsolution(Intersectmdbm)

    roots = map(x_sol, y_sol) do x,y
        [x,y]
    end
    
    f_vec(x_vec)::T = abs(disp(x_vec[1] + x_vec[2]*im))

    # select only those roots where the function is small
    fs = f_vec.(roots)
    fs = fs ./ mean(fs)

    inds1 = sortperm(fs)
    int1 = max(Int(round(length(inds1) * 0.2)), 4 * num_wavenumbers)
    inds1 = inds1[1:int1]

    w = max(0.2, abs(1.0 - std(fs)/3) + 10*tol)
    inds2 = findall(fs .< w)
    inds = intersect(inds1,inds2)

    roots = roots[inds]

    # refine the roots with optimisation
    roots = sort(roots, by = v -> v[2])

    # Used to create the NelderMead simplexer
    x_max = maximum(abs.(x_sol))
    y_max = maximum(abs.(y_sol))
    dx = T(0.5) * x_max / sqrt(length(roots))
    dy = T(0.5) * y_max / sqrt(length(roots))

    # Here we refine the roots
    roots2 = map(roots) do root

        inner_optimizer = NelderMead(
            initial_simplex =  EffectiveWaves.MySimplexer(dx,dy),
            parameters  = fixedparameters
        )

        res = optimize(f_vec, root, inner_optimizer, optimoptions)

        if res.minimum < T(100)*tol || (Optim.converged(res) && res.minimum < T(1000)*tol)
            res.minimizer
        else
            [zero(T),-one(T)]
        end
    end;
    deleteat!(roots2, findall(v-> v == [zero(T),-one(T)], roots2) )
    roots2 = reduce_kvecs(roots2, T(50) * tol * abs(kφ))

    k_effs = [sol[1] + sol[2]*im for sol in roots2]
    k_effs = sort(k_effs, by=imag)

    return k_effs
end
