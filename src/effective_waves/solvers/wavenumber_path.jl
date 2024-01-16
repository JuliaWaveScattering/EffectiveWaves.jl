# NOTE: PlanarAzimuthalSymmetry() does not include all possible wavenumbers
function wavenumbers_path(ω::T, micro::Microstructure{Dim};
        symmetry::AbstractSymmetry{Dim} = PlanarAzimuthalSymmetry{Dim}(),
        tol::T = 1e-5,
        mesh_points::Int = 3,
        mesh_size::Number = one(T) / T(2),
        num_wavenumbers::Int = 3,
        verbose::Bool = false,
        inner_optimizer = NelderMead(; parameters = NelderMeadparameters()),
        optimoptions::Optim.Options{T} = Optim.Options(g_tol = tol^T(3),x_tol=tol^T(2)),
        k_effs::Vector{Complex{T}} = Complex{T}[],
        # numberofparticles::T = Inf,
        kws...) where {T,Dim}

    medium = micro.medium

    # check parameters
        mesh_points = iseven(mesh_points) ? mesh_points + 1 : mesh_points

    # find at least one root to use as a scale for dk_x and dk_y
        eff_medium = effective_medium(micro)
        k0 = ω / eff_medium.c
        if isnan(k0) k0 = ω + T(0)*im end
        # abs(k0) can be used to non-dimensionlise k_vec
        # kscale = abs(k0)
        kscale = one(T)

        kφ = wavenumber_low_volumefraction(ω, micro;
            verbose = false
        )

        # guess initial mesh for lowest attenuating wavenumbers
        x_step = T(4) * mesh_size * max(abs(real(kφ)), real(k0), sqrt(eps(T)))
        xs = LinRange(-x_step,x_step,2 * mesh_points + 3)
        xs = [xs; [real(k0), real(kφ)]]

        # ys = [imag(kφ),imag(kφ)^2, sqrt(imag(kφ)), sqrt(eps(T))]
        ys = [imag(kφ), sqrt(eps(T))]
        ymin = minimum(ys);
        ymax = T(2) * mesh_size * maximum(ys)
        ys = LinRange(ymin, ymax, mesh_points)
        ys = [ys; [eps(T), imag(kφ)]]

        k_dim_vecs = [[x,y] for x in xs for y in ys]
        # k_vecs is non-dimensional
        k_vecs = k_dim_vecs ./ kscale

        low_tol = min(1e-4, sqrt(tol))

    # Generate some asymptotic wavenumbers
        k_asyms = asymptotic_monopole_wavenumbers(ω, micro;
            num_wavenumbers = num_wavenumbers + 2)
        k_asyms = k_asyms ./ kscale

    # Use asymptotic results to estimate maximum imaginary and real part.
        max_Imk = maximum([imag.(k_asyms); imag(kφ)])
        max_Rek = maximum([real.(k_asyms); real(kφ)])
        min_Rek = minimum([real.(k_asyms); real(kφ)])

        function constraint(vec::Vector{T})
            local cons = zero(T)
            if vec[1] > max_Rek
                cons = cons - one(T) + exp(T(10.0) * (vec[1] - max_Rek))
            end
            if vec[1] < min_Rek
                cons = cons - one(T) + exp(T(10.0) * (min_Rek - vec[1]))
            end
            if vec[2] > max_Imk
                cons = cons - one(T) + exp(T(10.0) * (vec[2] - max_Imk))
            end

            return cons
        end
        inds = findall(constraint.(k_vecs) .< one(T))
        k_vecs = k_vecs[inds]

    # The dispersion equation is given by: `dispersion([k1,k2]) = 0` where k_eff = k1 + im*k2.
        dispersion_dim = dispersion_equation(ω, micro, symmetry;
            tol = low_tol,
            kws...)

        dispersion(vec::Vector{T}) = constraint(vec) + dispersion_dim((vec[1] + vec[2]*im) .* kscale)

        k_vecs = [
            optimize(dispersion, kvec, #inner_optimizer,
            Optim.Options(x_tol=low_tol, g_tol = low_tol^3)
            ).minimizer
        for kvec in k_vecs]
        k_vecs = reduce_kvecs(k_vecs, low_tol/10)

        # add any specified keffs
        k_vecs = [
            [[real(kp),imag(kp)] for kp in k_asyms];
            k_vecs;
            [[real(kp),imag(kp)] ./ kscale for kp in k_effs]
        ]

        k_vecs = map(k_vecs) do k_vec
           res = optimize(dispersion, k_vec, inner_optimizer, optimoptions)
           if res.minimum < T(100)*tol || (Optim.converged(res) && res.minimum < T(10)*low_tol)
               res.minimizer
           else
               [zero(T),-one(T)]
           end
        end
        deleteat!(k_vecs, findall(v-> v == [zero(T),-one(T)], k_vecs) )
        k_vecs = reduce_kvecs(k_vecs, T(10)*tol)

        # Delete unphysical waves, including waves travelling backwards with almost no attenuation. This only is important in the limit of very low frequency or very weak scatterers.
        # deleteat!(k_vecs, findall([k_vec[2] < -sqrt(tol) for k_vec in k_vecs]))
        # delete wave travelling in wrong direction with small attenuation
        deleteat!(k_vecs, findall([-low_tol < abs(k_vec[2])/k_vec[1] < zero(T) for k_vec in k_vecs]))
        # deleteat!(k_vecs, find(k_vec[2] < tol && k_vec[1] < tol for k_vec in k_vecs))

    sort!(k_vecs, by= kv -> kv[2])

    # The best estimation of the mesh size is from asymptotic results

    dk = k_asyms[end] - k_asyms[end-1]

    dk_x = abs(real(dk)) * mesh_size
    dk_y = abs(imag(dk)) * mesh_size
    dk_xs = LinRange(-dk_x, dk_x, mesh_points)

    # Find two more roots, one with larger and smaller imaginary parts than the primary root kin
        two_roots = false # in case the roots are not on both branches of the tree, which is tested below
        kxs = dk_xs .+ k_vecs[1][1]
        ky =  k_vecs[1][2]
        invert = true

        # Find the first two roots that lead to the two branches of the root tree
        while !two_roots && (length(k_vecs) < num_wavenumbers) && ky <= max_Imk
            hits = [
                optimize(dispersion, [kx, ky], inner_optimizer, Optim.Options(x_tol = low_tol, g_tol = low_tol^3)).minimizer
            for kx in kxs]
            hits = reduce_kvecs(hits, low_tol)

            # Here we refine the hits
            hits = map(hits) do k_vec
                # res = optimize(dispersion, k_vec; g_tol = tol^2.0, f_tol = tol^4.0, x_tol=tol)
                res = optimize(dispersion, k_vec, inner_optimizer, optimoptions)
                if res.minimum < T(100)*tol || (Optim.converged(res) && res.minimum < T(10)*low_tol)
                    res.minimizer
                else
                    [zero(T),-one(T)]
                end
            end
            # deleteat!(hits, findall([-low_tol < abs(k_vec[2])/k_vec[1] < zero(T) for k_vec in hits]))
            deleteat!(hits, findall(v-> v == [zero(T),-one(T)], hits) )
            k_vecs = reduce_kvecs([hits;k_vecs], T(10)*tol)
            sort!(k_vecs, by= kv -> kv[2])

            if length(k_vecs) >= 2
                min1 = minimum(kv[1] for kv in k_vecs[2:end])
                max1 = maximum(kv[1] for kv in k_vecs[2:end])
                if min1*max1 < 0 || (max1-k_vecs[1][1])*(min1-k_vecs[1][1]) < 0
                    two_roots=true
                end
                dk_x += dk_x
                kxs = -k_vecs[2][1] .+ LinRange(-dk_x, dk_x, mesh_points)
                if invert
                    dk_y = k_vecs[2][2]/10
                    ky = k_vecs[2][2] - dk_y
                    invert = false
                elseif ky > 4*maximum(ks -> ks[2],k_vecs)
                    two_roots=true
                end
            elseif length(k_vecs) == 1
                dk_y += dk_y
            end
            ky = ky + dk_y
        end

    # Find more roots starting from the roots with the smaller imaginary part
        sort!(k_vecs, by = kv->kv[2])
        targets = k_vecs[2:end]
        dys = LinRange(0.5,1.0 + mesh_size, mesh_points)

        counter = 0
        while (num_wavenumbers > 2) && (counter < num_wavenumbers) && length(targets) > 0

            # procedure will be run at most once the desired number of wavenumbers
            counter = counter + 1

            # find roots following on from the smallest attenuating root
                sort!(targets, by = kv->kv[2])
                if verbose println("\n New target: $(targets[1])") end

            # find closest neighbour with smaller imaginary part
                i = findmin(
                    [(targets[1][2] != kvec[2]) ? norm(targets[1] - kvec) : Inf for kvec in k_vecs]
                )[2]
                # find at least one neighbour below, as the tree should grow
                j = (targets[1][2] > k_vecs[i][2]) ?
                    findmin(
                        [(targets[1][2] != kvec[2] && kvec[2] != k_vecs[i][2]) ? norm(targets[1] - kvec) : Inf for kvec in k_vecs]
                    )[2] :
                    findmin(
                        [(targets[1][2] > kvec[2] && kvec[2] != k_vecs[i][2]) ? norm(targets[1] - kvec) : Inf for kvec in k_vecs]
                    )[2]
                if j == Inf j = i end

                if verbose println("Closest neighbours: $(k_vecs[i]) and $(k_vecs[j])") end

            # create a mesh from this
                meshes = map([i,j]) do ind
                    tks = [k_vecs[ind] + dy.*(targets[1] - k_vecs[ind]) for dy in dys]
                    v = reverse(targets[1] - k_vecs[ind])
                    v[2] = -v[2]
                    hcat([
                        [tk + x.*v for x in LinRange(-mesh_size, mesh_size, mesh_points)]
                    for tk in tks]...)[:]
                end
                mesh = vcat(meshes...)
                deleteat!(mesh, findall([vec[2] < 0 || vec[2] > max_Imk for vec in mesh]))

            # search for roots from this mesh
                new_targets = map(mesh) do kin
                   optimize(dispersion, kin, inner_optimizer, Optim.Options(x_tol = low_tol, g_tol = low_tol^3)).minimizer
                end
                new_targets = reduce_kvecs(new_targets, low_tol/10)
                deleteat!(new_targets, findall(dispersion.(new_targets) .> low_tol))
                deleteat!(new_targets, findall([vec[2] < 0 || vec[2] > max_Imk for vec in new_targets]))

                # only keep targets which are not already in k_vecs
                new_targets = Vector{T}[
                        (findmin([norm(h - kvec) for kvec in k_vecs])[1] > low_tol) ? h : [zero(T),-one(T)]
                for h in new_targets]
                deleteat!(new_targets, findall(v-> v == [zero(T),-one(T)], new_targets))
                new_targets = reduce_kvecs(new_targets, low_tol)

                # Here we refine the new roots
                new_targets = map(new_targets) do k_vec
                    # res = optimize(dispersion, k_vec; g_tol = tol^2.0, f_tol = tol^4.0, x_tol=tol)
                    res = optimize(dispersion, k_vec, inner_optimizer, optimoptions)
                    if res.minimum < T(100)*tol || (Optim.converged(res) && res.minimum < T(10)*low_tol)
                        res.minimizer
                    else
                        [zero(T),-one(T)]
                    end
                end

                new_targets = [
                        (findmin([norm(h - kvec) for kvec in k_vecs])[1] > 10*tol) ? h : [zero(T),-one(T)]
                for h in new_targets]
                deleteat!(new_targets, findall(v-> v == [zero(T),-one(T)], new_targets))
                new_targets = reduce_kvecs(new_targets, T(10)*tol)

                if verbose println("New roots: $(new_targets)") end

                deleteat!(targets, 1)
                targets = Vector{T}[targets; new_targets]
                # group together wavenumbers which are closer than tol
                k_vecs = reduce_kvecs([new_targets; k_vecs], T(10)*tol)
                targets = filter(t -> t[2] <= max_Imk, targets)
        end

    # Finally delete unphysical waves, including waves travelling backwards with almost no attenuation. This only is important in the limit of very low frequency or very weak scatterers.
    # deleteat!(k_vecs, findall([k_vec[2] < -sqrt(tol) for k_vec in k_vecs]))
    # deleteat!(k_vecs, find(k_vec[2] < tol && k_vec[1] < tol for k_vec in k_vecs))
    deleteat!(k_vecs, findall([-low_tol < abs(k_vec[2])/k_vec[1] < zero(T) for k_vec in k_vecs]))

    # k_effs is dimensional
    k_effs = kscale .* [ k_vec[1] + k_vec[2]*im for k_vec in k_vecs]
    k_effs = sort(k_effs, by=imag)

    return k_effs
end
