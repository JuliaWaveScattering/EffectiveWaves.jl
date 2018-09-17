function wavenumbers_path(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; tol::T = 1e-6,
        hankel_order::Int = maximum_hankel_order(ω, medium, species; tol=tol),
        mesh_points::Int = 2, mesh_size::T = one(T),
        num_wavenumbers::Int = 3,
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; hankel_order = hankel_order),
        kws...) where T<:Number

    low_tol = 1e-4 # a tolerance used for a first pass with time_limit
    k = real(ω/medium.c)
    S = length(species)
    ho = hankel_order

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    function M_component(keff,j,l,m,n)
        (n==m ? 1.0:0.0)*(j==l ? 1.0:0.0) + 2.0pi*species[l].num_density*t_vecs[l][n+ho+1]*
            Nn(n-m,k*as[j,l],keff*as[j,l])/(k^2.0-keff^2.0)
    end

    # this matrix is needed to calculate the eigenvectors
    MM(keff::Complex{T}) = reshape(
        [M_component(keff,j,l,m,n) for m in -ho:ho, j = 1:S, n in -ho:ho, l = 1:S]
    , ((2ho+1)*S, (2ho+1)*S))

    # the constraint uses keff_vec[2] < -tol to better specify solutions where imag(k_effs)<0
    constraint(keff_vec::Array{T}) = ( (keff_vec[2] < -tol) ? one(T):zero(T))*(-1 + exp(-T(100.0)*keff_vec[2]))
    detMM2(keff_vec::Array{T}) =  constraint(keff_vec) + map(x -> real(x*conj(x)), det(MM(keff_vec[1]+im*keff_vec[2])))

    # find at least one root to use as a scale for dk_x and dk_y
        kφ = wavenumber_low_volfrac(ω, medium, species; verbose = false)
        eff_medium = effective_medium(medium, species)
        k0 = ω/eff_medium.c
        if isnan(k0) k0 = kφ end
        kin = [min(real(k0),abs(real(kφ))),abs(imag(kφ))]
        dx = kin[1]*mesh_size
        k_vecs = [[kin[1]+x,kin[2]] for x in linspace(-dx,dx,mesh_points+1)]
        push!(k_vecs,kin, [real(k0),zero(T)], [zero(T),zero(T)])
        k_vecs = [optimize(detMM2, kvec; x_tol= low_tol).minimizer for kvec in k_vecs]
        k_vecs = reduce_kvecs(k_vecs, low_tol/20)

    sort!(k_vecs, by= kv -> kv[2])
    dk_x = abs(k_vecs[1][1]) * mesh_size
    dk_y = abs(k_vecs[1][2]) * mesh_size
    dk_xs = linspace(-dk_x, dk_x, mesh_points)

    # Find two more roots, one with larger and smaller imaginary parts than the primary root kin
        two_roots=false
        kxs = dk_xs .+ k_vecs[1][1]
        ky =  k_vecs[1][2]
        invert = true

        while !two_roots && (length(k_vecs) < num_wavenumbers)
            ky = ky + dk_y
            hits = [
                optimize(detMM2, [kx, ky]; x_tol = low_tol).minimizer
            for kx in kxs]
            deleteat!(hits, find(detMM2.(hits) .> low_tol))
            k_vecs = reduce_kvecs([hits;k_vecs], low_tol/20)
            sort!(k_vecs, by= kv -> kv[2])

            if length(k_vecs) > 2
                sort!(k_vecs, by= kv -> kv[2])
                min1 = minimum(kv[1] for kv in k_vecs[2:end])
                max1 = maximum(kv[1] for kv in k_vecs[2:end])
                if min1*max1 < 0 || (max1-k_vecs[1][1])*(min1-k_vecs[1][1]) < 0
                    two_roots=true
                end
                dk_x += dk_x
                kxs = -k_vecs[2][1] + linspace(-dk_x, dk_x, mesh_points)
            elseif length(k_vecs) == 1
                dk_y += dk_y
            elseif invert
                sort!(k_vecs, by= kv -> kv[2])
                kxs = dk_xs .- k_vecs[2][1]
                dk_y = k_vecs[2][2]/10
                ky = k_vecs[2][2] - dk_y
                invert = false
            elseif ky > 4*maximum(ks -> ks[2],k_vecs)
                two_roots=true
            else
                dk_x += dk_x
                kxs = -k_vecs[2][1] + linspace(-dk_x, dk_x, mesh_points)
            end
        end

    # Find roots following on from the two roots
        sort!(k_vecs, by = kv->kv[2])
        targets = k_vecs[2:end]
        dys = linspace(0.5*mesh_size,1.0+mesh_size,mesh_points)
        while (length(k_vecs) < num_wavenumbers) && length(targets) > 0
            # find roots following on from the smallest attenuating root
                sort!(targets, by = kv->kv[2])

            # find closest neighbour with smaller imaginary part
                i = findmin((targets[1][2] != kvec[2])? norm(targets[1] - kvec) : Inf for kvec in k_vecs)[2]
                j = findmin((targets[1][2] != kvec[2] && kvec[2] != k_vecs[i][2])? norm(targets[1] - kvec) : Inf for kvec in k_vecs)[2]

            # create a mesh from this
                meshes = map([i,j]) do ind
                    tks = [k_vecs[ind] + dy.*(targets[1] - k_vecs[ind]) for dy in dys]
                    v = reverse(targets[1] - k_vecs[ind])
                    v[2] = -v[2]
                    hcat([
                        [tk + x.*v for x in linspace(-mesh_size, mesh_size, mesh_points)]
                    for tk in tks]...)[:]
                end
                mesh = vcat(meshes...)
                deleteat!(mesh, find(vec[2] < 0 for vec in mesh))

            # search for roots from this mesh
                hits = map(mesh) do kin
                   optimize(detMM2, kin; x_tol = low_tol, f_tol = low_tol/10).minimizer
                end
                hits = reduce_kvecs(hits, low_tol/20)
                deleteat!(hits, find(detMM2.(hits) .> low_tol))

            # find new targets
                new_targets = [
                        (findmin([norm(h - kvec) for kvec in k_vecs])[1] > low_tol)? h : [zero(T),-one(T)]
                    for h in hits]

                deleteat!(targets, 1)
                targets = [targets; new_targets]
                deleteat!(targets, find(vec[2] < 0 for vec in targets))
                k_vecs = reduce_kvecs([hits;k_vecs], low_tol/20)
        end

    # Here we refine the effective wavenumbers
    k_vecs = map(k_vecs) do k_vec    # Here we refine the effective wavenumbers
        res = optimize(detMM2, k_vec; g_tol = tol^2.0, f_tol = tol^4.0, x_tol=tol)
        if res.minimum > T(20)*tol
            [zero(T),-one(T)]
        else
            res.minimizer
        end
    end

    # group together wavenumbers which are closer than tol
    k_vecs = reduce_kvecs(k_vecs,tol)

    # Finally delete unphysical waves, including waves travelling backwards with almost no attenuation. This only is important in the limit of very low frequency or very weak scatterers.
    deleteat!(k_vecs, find(k_vec[2] < -T(10)*tol for k_vec in k_vecs))
    deleteat!(k_vecs, find(k_vec[2] < tol && k_vec[1] < tol for k_vec in k_vecs))

    k_effs = [ k_vec[1] + k_vec[2]*im for k_vec in k_vecs]
    k_effs = sort(k_effs, by=imag)

    return k_effs
end
