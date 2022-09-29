
# convert the PlanarAzimuthalSymmetry eigenvector to the PlanarSymmetry eigenvector
function convert_eigenvector_basis(medium::Acoustic{T,3},sym::PlanarAzimuthalSymmetry{3},eigvecs::Array{Complex{T}}) where T
    basis_order = size(eigvecs,1) - 1
    S = size(eigvecs,2)
    P = size(eigvecs,3)

    v = [
        (m == 0) ? eigvecs[l+1,s,p] : zero(Complex{T})
    for p = 1:P for s = 1:S for l = 0:basis_order for m = -l:l]

    return reshape(v,basisorder_to_basislength(Acoustic{T,3},basis_order),S,P)
end

function eigensystem(ω::T, medium::Acoustic{T,2}, micro::ParticulateMicrostructure{2}, ::AbstractPlanarSymmetry;
        basis_order::Int = 2,
        numberofparticles::Number = Inf,
        kws...) where {T<:AbstractFloat}

    scale_number_density = one(T)
    if numberofparticles > 1
        scale_number_density = one(T) - one(T) / numberofparticles
    end

    k = ω / medium.c
    sps = micro.species
    S = length(sps)
    ho = basis_order

    t_matrices = get_t_matrices(medium, sps, ω, ho)

    len = (2ho+1) * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    as = [
        s1.exclusion_distance * outer_radius(s1) + s2.exclusion_distance * outer_radius(s2)
    for s1 in sps, s2 in sps]

    if length(micro.paircorrelations[1].r) > 1
        pair_rs, hks, gs = precalculate_pair_correlations(micro, k, ho)
    end

    function M_component(keff,Ns,j,l,m,n)
        (n == m ? 1.0 : 0.0)*(j == l ? 1.0 : 0.0) + 2.0pi * scale_number_density * number_density(sps[l]) * t_matrices[l][m+ho+1,m+ho+1] * Ns[n-m + 2ho+1,j,l] / (keff^2.0 - k^2.0)
    end

    # this matrix is needed to calculate the eigenvectors
    function MM(keff::Complex{T})
        Ns = [
            kernelN2D(m,k*as[s1,s2],keff*as[s1,s2])
        for m = (-2ho):2ho, s1 = 1:S, s2 = 1:S]

        # For a pair correlation which is not hole correction need to add a finite integral
        if length(micro.paircorrelations[1].r) > 1
            Ns = Ns + kernelW2D(k, keff, pair_rs, gs, hks, basis_order)
        end

        ind2 = 1
        for l = 1:S for n = -ho:ho
            ind1 = 1
            for s = 1:S for m = -ho:ho
                MM_mat[ind1,ind2] = M_component(keff,Ns,s,l,m,n)
                ind1 += 1
        end end
            ind2 += 1
        end end
        return MM_mat
    end

    return MM
end

function eigensystem(ω::T, medium::Acoustic{T,3}, micro::ParticulateMicrostructure{3}, ::AbstractPlanarSymmetry;
        basis_order::Int = 2,
        direction_eff::Union{AbstractVector{T},AbstractVector{Complex{T}}} = [0.0,0.0,1.0],
        numberofparticles::Number = Inf,
        kws...) where {T<:AbstractFloat}

    scale_number_density = one(T)
    if numberofparticles > 1
        scale_number_density = one(T) - one(T) / numberofparticles
    end

    species = micro.species

    k = ω/medium.c
    S = length(species)
    ho = basis_order
    len = basisorder_to_basislength(Acoustic{T,3},ho) * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    rθφp = cartesian_to_radial_coordinates(direction_eff)
    # below we take the real part because our associated Legendre functions are only implemented for real arguments.
    if norm(real.(rθφp) - rθφp) / norm(rθφp) < (eps(T))^(1/3)
        rθφp = real.(rθφp)
    end

    Ys = spherical_harmonics(2ho, rθφp[2], rθφp[3]);
    lm_to_n = lm_to_spherical_harmonic_index

    t_matrices = get_t_matrices(medium, species, ω, ho)
    t_diags = diag.(t_matrices)
    baselen(order::Int) = basisorder_to_basislength(Acoustic{T,3},order)

    as = [
        s1.exclusion_distance * outer_radius(s1) + s2.exclusion_distance * outer_radius(s2)
    for s1 in species, s2 in species]

    if length(micro.paircorrelations[1].r) > 1
        pair_rs, hks, gs = precalculate_pair_correlations(micro, k, ho)
    end

    function M_component(keff::Complex{T},Ns::Array{Complex{T}},l::Int,m::Int,s1::Int,dl::Int,dm::Int,s2::Int)::Complex{T}
        (m == dm && l == dl && s1 == s2 ? one(Complex{T}) : zero(Complex{T})) +
        4pi * as[s1,s2] * scale_number_density * number_density(species[s2]) * t_diags[s1][baselen(l)] *
        sum(
            Complex{T}(im)^(-l1) * Ys[lm_to_n(l1,dm-m)] * Ns[l1+1,s1,s2] *
            gaunt_coefficient(dl,dm,l,m,l1,dm-m)
        for l1 in max(abs(dm-m),abs(dl-l)):(dl+l)) / (keff^2.0 - k^2.0)
    end

    function MM(keff::Complex{T})::Matrix{Complex{T}}
        Ns = [kernelN3D(l,k*as[s1,s2],keff*as[s1,s2]) for l = 0:2ho, s1 = 1:S, s2 = 1:S]

        # For a pair correlation which is not hole correction need to add a finite integral
        if length(micro.paircorrelations[1].r) > 1
            Ns = Ns + kernelW3D(k, keff, pair_rs, gs, hks, basis_order)
        end

        ind2 = 1
        for s2 = 1:S, dl = 0:ho for dm = -dl:dl
            ind1 = 1
            for s1 = 1:S, l = 0:ho for m = -l:l
                MM_mat[ind1, ind2] = M_component(keff,Ns,l,m,s1,dl,dm,s2)
                ind1 += 1
            end end
            ind2 += 1
        end end
        return MM_mat
    end

    return MM
end

function eigensystem(ω::T, medium::Acoustic{T,3}, micro::ParticulateMicrostructure{3}, ::PlanarAzimuthalSymmetry{3};
        basis_order::Int = 2,
        numberofparticles::Number = Inf,
        kws...) where {T<:AbstractFloat}

    scale_number_density = one(T)
    if numberofparticles > 1
        scale_number_density = one(T) - one(T) / numberofparticles
    end

    species = micro.species

    k = ω/medium.c
    S = length(species)
    ho = basis_order
    len = (ho+1) * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    t_matrices = get_t_matrices(medium, species, ω, ho)
    t_diags = diag.(t_matrices)
    baselen(order::Int) = basisorder_to_basislength(Acoustic{T,3},order)

    as = [
        s1.exclusion_distance * outer_radius(s1) + s2.exclusion_distance * outer_radius(s2)
    for s1 in species, s2 in species]

    # Pre calculations for pair correlation
    # Am going to assume the discrete pair correlation is sampled on the same mesh for every specie. Otherwise the code will be too inefficient.
    if length(micro.paircorrelations[1].r) > 1
        pair_rs, hks, gs = precalculate_pair_correlations(micro, k, ho)
    end

    function M_component(keff::Complex{T}, Ns::Array{Complex{T}},l,s1,dl,s2)::Complex{T}
        (l == dl && s1 == s2 ? one(Complex{T}) : zero(Complex{T})) +
        4pi * as[s1,s2] * scale_number_density * number_density(species[s2]) * t_diags[s1][baselen(l)] *
        sum(
            Complex{T}(im)^(-l1) * sqrt((2*l1+1)/(4pi) ) * Ns[l1+1,s1,s2] *
            gaunt_coefficient(dl,0,l,0,l1,0)
        for l1 in abs(dl-l):(dl+l)) ./ (keff^2.0 - k^2.0)
    end

    function MM(keff::Complex{T})::Matrix{Complex{T}}

        Ns = [
            kernelN3D(l,k*as[s1,s2],keff*as[s1,s2])
        for l = 0:2ho, s1 = 1:S, s2 = 1:S]

        # For a pair correlation which is not hole correction need to add a finite integral
        if length(micro.paircorrelations[1].r) > 1
            Ns = Ns + kernelW3D(k, keff, pair_rs, gs, hks, basis_order)
        end

        ind2 = 1
        for s2 = 1:S, dl = 0:ho
            ind1 = 1
            for s1 = 1:S, l = 0:ho
                MM_mat[ind1, ind2] = M_component(keff,Ns,l,s1,dl,s2)
                ind1 += 1
            end
            ind2 += 1
        end
        return MM_mat
    end

    return MM
end
