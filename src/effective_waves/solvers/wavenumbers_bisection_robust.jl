wavenumbers_bisection_robust(ω::Number, medium::PhysicalMedium, specie::Specie; kws...) = wavenumbers_bisection_robust(ω, Microstructure(medium,[specie]); kws...)
wavenumbers_bisection_robust(ω::Number, medium::PhysicalMedium, species::Species; kws...) = wavenumbers_bisection_robust(ω, Microstructure(medium,species); kws...)

# NOTE: PlanarAzimuthalSymmetry() does not included all possible wavenumbers
function wavenumbers_bisection_robust(ω::T, micro::Microstructure{Dim};
        symmetry::AbstractSymmetry{Dim} = PlanarAzimuthalSymmetry{Dim}(),
        basis_order::Int = 2,
        tol::T = 1e-5,
        distance_factor::T = 1.0,
        num_wavenumbers = 1,
        box_k::Vector{Vector{T}} = box_keff(ω, micro; num_wavenumbers = num_wavenumbers, basis_order = basis_order),
        bisection_mesh_points::Int = 2 * Int(round(- log(tol))) + 2*num_wavenumbers,
        bisection_iteration::Int = 2,
        fixedparameters::Optim.FixedParameters = NelderMeadparameters(),
        optimoptions::Optim.Options{T} = Optim.Options(
            iterations =  Int(round(-log(tol))) * 200,
            g_tol = tol^T(3), x_abstol=tol^T(2)),
        kws...) where {T,Dim}

    kφ = wavenumber_low_volumefraction(ω, micro;
        verbose = false, basis_order = basis_order,
    )
    disp = dispersion_complex(ω, micro, symmetry; basis_order = basis_order, kws...)

    freal(x, y)::T = real(disp(x + y*im))
    fimag(x, y)::T = imag(disp(x + y*im))

    x = LinRange(box_k[1][1],box_k[1][2],bisection_mesh_points)
    x = [x; -real(kφ):real(kφ):(real(kφ))]
    x = sort(x)
    y = LinRange(box_k[2][1],box_k[2][2],bisection_mesh_points)
    y = [y; -imag(kφ):imag(kφ):(imag(kφ))]
    y = sort(y)

    # axes=[range ω, real part of kz, imag part of kz]
    # axes = [1:3,0:4,-4:4]

    # Doesn't work as real part has tiny loop around the principal root.
        # #Intersection of both real and imag surface of dispersion equation
        #     axes = [Axis(ω_mesh,"ω"), Axis(x,"x"),Axis(y,"y")]
        #     f_combined(x...) = (freal(x...), fimag(x...))
        #     Intersectmdbm = MDBM_Problem(f_combined,axes)

    axes_xy = [Axis(x,"x"),Axis(y,"y")];
    
    Intersectmdbm = MDBM_Problem(fimag,axes_xy);
    solve!(Intersectmdbm,bisection_iteration);
    x_sol_imag, y_sol_imag = getinterpolatedsolution(Intersectmdbm);
    
    Intersectmdbm = MDBM_Problem(freal,axes_xy)

    solve!(Intersectmdbm,bisection_iteration)
    x_sol_real, y_sol_real = getinterpolatedsolution(Intersectmdbm)

    x_sol = vcat(x_sol_imag, x_sol_real)
    y_sol = vcat(y_sol_imag, y_sol_real)

    roots = map(x_sol, y_sol) do x,y
        [x,y]
    end

    f_vec(x_vec)::T = abs(disp(x_vec[1] + x_vec[2]*im))

    # select only 10% of the smallest roots
    fs = f_vec.(roots)
    fs = fs ./ mean(fs)

    inds1 = sortperm(fs)
    int1 = max(Int(round(length(inds1) * 0.1)), 4 * num_wavenumbers)
    inds1 = inds1[1:int1]

    # select only those roots where the function is small
    # println("abs(1.0 - std(fs)/2): ", abs(1.0 - std(fs)/2))

    w = max(0.2, abs(1.0 - std(fs)) + 10*tol)
    inds2 = findall(fs .< w)
    
    inds = intersect(inds1,inds2)
    roots = roots[inds]

    # expected distance between roots:
    k_asyms = asymptotic_monopole_wavenumbers(ω, micro;
    num_wavenumbers = 2)
    k_asyms_vecs = [[real(k), imag(k)] for k in k_asyms]

    dist = if k_asyms_vecs[1] ∈ box_k || k_asyms_vecs[2] ∈ box_k
        abs(k_asyms[3] - k_asyms[2]) * 0.6
    else
        Inf
    end
    
    # the most important distance is the between the main roots and its negative counterpart, which is 2*abs(kφ). This distance also tends to be smaller then the dist above.
    dist = distance_factor * min(1.5*abs(kφ), dist)

    # find clusters of roots that are close together, and replace them with the smallest roots.
    all_inds = collect(eachindex(roots))
    roots_vec = map(roots) do r
        ind_ins = findall([norm(v - r) < dist for v in roots[all_inds]])
        inds = all_inds[ind_ins]
        deleteat!(all_inds,ind_ins)

        if isempty(inds)
            [[zero(T),-one(T)]]
        else
            finds = sortperm(f_vec.(roots[inds]))
            roots[inds[
                finds[1:min(2, length(finds))]
            ]]
        end
    end
    deleteat!(roots_vec, findall(v-> v == [[zero(T),-one(T)]], roots_vec) )

    roots = vcat(roots_vec...)
    roots = reduce_kvecs(roots, tol * abs(kφ))
    roots = [roots; [[real(kφ), imag(kφ)]]]

    # # Used to create the NelderMead simplexer
    # x_max = maximum(abs.(x_sol))
    # y_max = maximum(abs.(y_sol))
    # dx = T(0.5) * x_max / sqrt(length(roots))
    # dy = T(0.5) * y_max / sqrt(length(roots))

    # Here we refine the roots
    roots2 = map(roots) do root
        inner_optimizer = NelderMead(
            # initial_simplex =  EffectiveWaves.MySimplexer(dx,dy),
            parameters  = fixedparameters
        )

        res = optimize(f_vec, root, inner_optimizer, optimoptions)

        if res.minimum < sqrt(tol) || (Optim.converged(res) && res.minimum < T(10)* sqrt(tol))
            res.minimizer
        else
            [zero(T),-one(T)]
        end
    end;
    deleteat!(roots2, findall(v-> v == [zero(T),-one(T)], roots2) )
    roots2 = reduce_kvecs(roots2, T(50) * tol * abs(kφ))

    low_tol =  T(50) * tol * abs(kφ)
    inds = findall([-low_tol < abs(k_vec[2])/k_vec[1] < zero(T) for k_vec in roots2])
    deleteat!(roots2, inds)

    k_effs = [sol[1] + sol[2]*im for sol in roots2]
    k_effs = sort(k_effs, by=imag)

    return k_effs
end
