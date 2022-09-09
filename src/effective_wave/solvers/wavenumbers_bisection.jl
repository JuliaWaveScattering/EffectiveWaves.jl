# NOTE: PlanarAzimuthalSymmetry() does not included all possible wavenumbers
function wavenumbers_bisection(ω::T, medium::PhysicalMedium{Dim}, species::Species{Dim};
        symmetry::AbstractSymmetry{Dim} = PlanarAzimuthalSymmetry{Dim}(),
        tol::T = 1e-5,
        num_wavenumbers = 3,
        box_k::Vector{Vector{T}} = box_keff(ω, medium, species; tol = tol),
        bisection_mesh_points::Int = Int(round(2 - log(tol))) + 25 * num_wavenumbers,
        fixedparameters::Optim.FixedParameters = NelderMeadparameters(),
        optimoptions::Optim.Options{T} = Optim.Options(
            g_tol = tol^T(3), x_tol=tol^T(2),
            iterations = 800 - Int(round(log(tol)*40))
        ),
        kws...) where {T,Dim}

    ko = real(ω / medium.c)

    disp = dispersion_complex(ω, medium, species, symmetry; kws...)
    # disp = dispersion_complex(ω, medium, species, symmetry; basis_order = basis_order)

    function f(x, y)
        z = disp(x + y*im)
        [real(z),imag(z)]
    end

    x = LinRange(box_k[1][1],box_k[1][2],bisection_mesh_points);
    y = LinRange(box_k[2][1],box_k[2][2], Int(round(bisection_mesh_points/2)));

    mymdbm = MDBM_Problem(f,[Axis(x,"x"),Axis(y,"y")]);

    iteration = 1 #number of refinements (resolution doubling)
    solve!(mymdbm,iteration; interpolationorder = 0);

    x_sol,y_sol=getinterpolatedsolution(mymdbm);
    sols = map(x_sol,y_sol) do x,y
        [x,y]
    end
    sols = reduce_kvecs(sols, ko * sqrt(tol) * T(10))

    # using Plots
    # x_eval,y_eval=getevaluatedpoints(mymdbm)
    # scatter(x_eval[1:7:end],y_eval[1:7:end], label="eval", markersize = 0.0, markerstrokewidth = 0.0)
    # scatter!(x_sol,y_sol, label="sol", markersize = 1.0, markerstrokewidth = 2.0)

    # keep only those with the smallest attenuation
    sols = sort(sols, by = vec -> vec[2])
    sols = sols[1:min(length(sols),10*num_wavenumbers)]

    # Add some guess from asymptotics
    kφ = wavenumber_low_volumefraction(ω, medium, species)

    k_real = max(abs(real(kφ)), real(ko), sqrt(eps(T)))
    k_imag = [imag(kφ), sqrt(eps(T))]

    k_vecs = [[x,y] for x in LinRange(-k_real,k_real,10) for y in k_imag]
    sols = [sols; k_vecs]

    function f_vec(vec::Vector{T})
        abs(disp(vec[1] + vec[2]*im))
    end

    dx = T(0.5) * maximum(abs.(x_sol)) / sqrt(length(sols))
    dy = T(0.5) * maximum(abs.(y_sol)) / sqrt(length(sols))

    # Here we refine the roots
    sols = map(sols) do vec
        inner_optimizer = NelderMead(
            initial_simplex =  MySimplexer(dx,dy),
            parameters  = fixedparameters
        )

        res = optimize(f_vec, vec, inner_optimizer, optimoptions)
        if res.minimum < T(10)*tol || (Optim.converged(res) && res.minimum < T(100)*tol)
            res.minimizer
        else
            [zero(T),-one(T)]
        end
    end
    sols = deleteat!(sols, findall(v-> v == [zero(T),-one(T)], sols))
    sols = reduce_kvecs(sols, T(10) * tol * ko)

    k_effs = [sol[1] + sol[2]*im for sol in sols]
    k_effs = sort(k_effs, by=imag)

    return k_effs
end
