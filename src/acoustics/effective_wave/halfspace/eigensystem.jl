function eigensystem(ω::T, medium::PhysicalMedium{T,2}, species::Species{T,2}, ::AbstractPlanarSymmetry;
        basis_order::Int = 2,
        kws...) where {T<:AbstractFloat}

    k = ω / medium.c
    sps = species
    S = length(sps)
    ho = basis_order

    t_matrices = get_t_matrices(medium, sps, ω, ho)

    len = (2ho+1) * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    as = [
        s1.exclusion_distance * outer_radius(s1) + s2.exclusion_distance * outer_radius(s2)
    for s1 in sps, s2 in sps]

    function M_component(keff,j,l,m,n)
        (n == m ? 1.0 : 0.0)*(j == l ? 1.0 : 0.0) - 2.0pi * number_density(sps[l]) * t_matrices[l][m+ho+1,m+ho+1] * kernelN2D(n-m,k*as[j,l],keff*as[j,l]) / (k^2.0-keff^2.0)
    end

    # this matrix is needed to calculate the eigenvectors
    function MM(keff::Complex{T})
        ind2 = 1
        for l = 1:S for n = -ho:ho
            ind1 = 1
            for s = 1:S for m = -ho:ho
                MM_mat[ind1,ind2] = M_component(keff,s,l,m,n)
                ind1 += 1
        end end
            ind2 += 1
        end end
        return MM_mat
    end

    return MM
end

function eigensystem(ω::T, medium::PhysicalMedium{T,3}, species::Species{T,3}, ::AbstractPlanarSymmetry;
        basis_order::Int = 2,
        θp::Complex{T} = zero(Complex{T}),
        φp::Complex{T} = zero(Complex{T}),
        kws...) where {T<:AbstractFloat}

    k = real(ω/medium.c)
    S = length(species)
    ho = basis_order
    len = (ho+1)^2 * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    Ys = spherical_harmonics(L1, θp, φp);
    lm_to_n = lm_to_spherical_harmonic_index

    t_matrices = get_t_matrices(medium, sps, ω, ho)

    as = [ s1.exclusion_distance * outer_radius(s1) + s2.exclusion_distance * outer_radius(s2) for s1 in species, s2 in species]
    function M_component(keff::Complex{T},Ns::Array{Complex{T}},l::Int,m::Int,s1::Int,dl::Int,dm::Int,s2::Int)::Complex{T}
        (m == dm && l == dl && s1 == s2 ? one(Complex{T}) : zero(Complex{T})) +
        4pi * as[s1,s2] * number_density(species[s2]) * t_matrices[s1][l+ho+1,l+ho+1] *
        sum(
            Complex{T}(im)^(-l1) * Ys[lm_to_n(l1,dm-m)] * Ns[l1+1,s1,s2] *
            gaunt_coefficients(dl,dm,l,m,l1,dm-m)
        for l1 in max(abs(dm-m),abs(dl-l)):(dl+l)) / (keff^2.0 - k^2.0)
    end

    function MM(keff::Complex{T})::Matrix{Complex{T}}
        Ns = [kernelN3D(l,k*as[s1,s2],keff*as[s1,s2]) for l = 0:2ho, s1 = 1:S, s2 = 1:S]
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

function eigensystem(ω::T, medium::Acoustic{T,3}, species::Species{T}, ::PlanarAzimuthalSymmetry;
        basis_order::Int = 2,
        kws...) where {T<:AbstractFloat}



    k = real(ω/medium.c)
    S = length(species)
    ho = basis_order
    len = (ho+1) * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    t_matrices = get_t_matrices(medium, sps, ω, ho)

    as = [ s1.exclusion_distance * outer_radius(s1) + s2.exclusion_distance * outer_radius(s2) for s1 in species, s2 in species]
    function M_component(keff::Complex{T},Ns::Array{Complex{T}},l,s1,dl,s2)::Complex{T}
        (l == dl && s1 == s2 ? one(Complex{T}) : zero(Complex{T})) +
        4pi * as[s1,s2] * number_density(species[s2]) * t_matrices[s1][l+ho+1,l+ho+1] *
        sum(
            Complex{T}(im)^(-l1) * sqrt((2*l1+1)/(4pi) ) * Ns[l1+1,s1,s2] *
            gaunt_coefficients(dl,0,l,0,l1,0)
        for l1 in abs(dl-l):(dl+l)) / (keff^2.0 - k^2.0)
    end

    function MM(keff::Complex{T})::Matrix{Complex{T}}
        Ns = [kernelN3D(l,k*as[s1,s2],keff*as[s1,s2]) for l = 0:2ho, s1 = 1:S, s2 = 1:S]
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
