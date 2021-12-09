# NOTE: PlanarAzimuthalSymmetry() does not included all possible wavenumbers
function wavenumbers_bisection_robust(ω::T, medium::PhysicalMedium{Dim}, species::Species{T,Dim};
        symmetry::AbstractSymmetry{Dim} = PlanarAzimuthalSymmetry{Dim}(),
        tol::T = 1e-5,
        num_wavenumbers = 3,
        box_k::Vector{Vector{T}} = box_keff(ω, medium, species; tol = tol),
        bisection_mesh_points::Int = Int(round(2 - log(tol))) + 3 * num_wavenumbers,
        bisection_iteration::Int = Int(round(-log(tol) / 3)),
        fixedparameters::Optim.FixedParameters = NelderMeadparameters(),
        optimoptions::Optim.Options{T} = Optim.Options(
            g_tol = tol^T(3), x_tol=tol^T(2),
            iterations = 1000 + Int(round(log(tol)*20))
        ),
        kws...) where {T,Dim}

    # exact low-frequency effective wavenumber is known exactly
    ko = real(ω / effective_medium(medium, species).c)

    disp = dispersion_complex(ω, medium, species, symmetry; kws...)
    # disp = dispersion_complex(ω, medium, species, symmetry; basis_order = basis_order)

    freal(x, y)::T = real(disp(x + y*im))
    fimag(x, y)::T = imag(disp(x + y*im))

    x = LinRange(box_k[1][1],box_k[1][2],bisection_mesh_points);
    y = LinRange(box_k[2][1],box_k[2][2], Int(round(bisection_mesh_points/2)));

    prob_real = MDBM_Problem(freal,[Axis(x,"x"),Axis(y,"y")]);
    prob_imag = MDBM_Problem(fimag,[Axis(x,"x"),Axis(y,"y")]);

    solve!(prob_real,bisection_iteration);
    solve!(prob_imag,bisection_iteration);

    # x_eval,y_eval=getevaluatedpoints(prob_real)
    xr_sol,yr_sol=getinterpolatedsolution(prob_real);
    xi_sol,yi_sol=getinterpolatedsolution(prob_imag)

    x_max = maximum(abs.([xr_sol; xi_sol]))
    y_max = maximum(abs.([yr_sol; yi_sol]))

    x_normalise = sqrt(length(xr_sol)) / x_max
    y_normalise = sqrt(length(yr_sol)) / y_max

    vecs1 = map(xr_sol,yr_sol) do x,y
        round.([x * x_normalise, y * y_normalise])
    end
    vecs2 = map(xi_sol,yi_sol) do x,y
        round.([x * x_normalise, y * y_normalise])
    end

    root_vecs = collect(intersect(Set(vecs1),Set(vecs2)))

    inds1 = [
        findall([v == root for v in vecs1])
    for root in root_vecs]

    inds2 = [
        findall([v == root for v in vecs2])
    for root in root_vecs]

    roots = map(eachindex(root_vecs)) do i
        x = mean([xr_sol[inds1[i]];xi_sol[inds2[i]]])
        y = mean([yr_sol[inds1[i]];yi_sol[inds2[i]]])
        [x, y]
    end
    roots = reduce_kvecs(roots, T(50) * tol * ko)

    # If the mesh is not fine enough,
    roots = sort(roots, by = v -> v[2])
    roots = roots[1:min(length(roots),50*num_wavenumbers)]

    # using Plots
    #points where the function foo was evaluated
    # x_eval,y_eval=getevaluatedpoints(prob_real)
    # scatter(x_eval,y_eval, label="real", markersize = 1.0, markerstrokewidth = 0.0)

    # scatter(xr_sol,yr_sol, label="real", markersize = 1.0, markerstrokewidth = 0.0)
    # scatter!(xi_sol,yi_sol, label="imag", markersize = 1.0, markerstrokewidth = 0.0)
    # scatter!([r[1] for r in roots], [r[2] for r in roots], markersize = 2.0)

    f_vec(x_vec)::T = abs(disp(x_vec[1] + x_vec[2]*im))

    # Used to create the NelderMead simplexer
    dx = T(0.5) * x_max / sqrt(length(roots))
    dy = T(0.5) * y_max / sqrt(length(roots))

    # Here we refine the roots
    roots2 = map(roots) do root

        inner_optimizer = NelderMead(
            initial_simplex =  MySimplexer(dx,dy),
            parameters  = fixedparameters
        )

        res = optimize(f_vec, root, inner_optimizer, optimoptions)

        if res.minimum < T(100)*tol || (Optim.converged(res) && res.minimum < T(1000)*tol)
            res.minimizer
        else
            [zero(T),-one(T)]
        end
    end
    roots2 = reduce_kvecs(roots2, T(50) * tol * ko)

    k_effs = [sol[1] + sol[2]*im for sol in roots2]
    k_effs = sort(k_effs, by=imag)

    return k_effs
end
