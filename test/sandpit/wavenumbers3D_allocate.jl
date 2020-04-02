
function wavematrix3D_allocate(ω::T, medium::Medium{T}, species::Species{T};
        dim = 3, tol::T = 1e-4,
        basis_order::Int = maximum_basis_order(ω, medium, species; tol=tol),
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; basis_order = basis_order, dim = dim),
        kws...) where T<:Number

    k = real(ω/medium.c)
    S = length(species)
    ho = basis_order
    len = (ho+1)^4 * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    cs = reshape(
        [
            gaunt_coefficients(l1,m1,l2,m2,l3,m3)
        for l1 = 0:ho for m1 = -l1:l1 for l2 = 0:ho for m2 = -l2:l2 for l3 = 0:ho for m3 = -l3:l3]
    , ((ho+1)^2,(ho+1)^2,(ho+1)^2));

    # NOTE that cs[n3,n2,n1] == gaunt_coefficients(l1,m1,l2,m2,l3,m3)
    ls, ms = spherical_harmonics_indices(ho);
    n_max = length(ls)

    function M_component2(keff,Ns,n,dn,n1,n2,s1,s2)::Complex{T}
        (n == dn && n1 == n2 && s1 == s2 ? 1.0 : 0.0) +
        as[s1,s2] * species[s2].num_density * t_vecs[s1][ls[n]+ho+1] *
        sum(cs[:,n,dn] .* cs[:,n1,n2] .* Ns[ls[:].+1,s1,s2]) / (keff^2.0 - k^2.0)
    end
    function MM2(keff::Complex{T})
        Ns = [kernelN2D(l3,k*as[s1,s2],keff*as[s1,s2]; dim = dim) for l3 = 0:ho, s1 = 1:S, s2 = 1:S]
        ind2 = 1
        for s2 = 1:S for dn = 1:n_max for n1 = 1:n_max
            ind1 = 1
            for s1 = 1:S for n = 1:n_max for n2 = 1:n_max
                MM_mat[ind1, ind2] = M_component2(keff,Ns,n,dn,n1,n2,s1,s2)
                ind1 += 1
            end end end
            ind2 += 1
        end end end
        return MM_mat
    end
    # keff = 1.0 + rand()*im
    # s1 = s2 = 1

    # norm([
    #     M_component2(keff,Ns,n,dn,n1,n2,s1,s2) - M_component(keff,Ns,ls[n],ms[n],ls[n2],ms[n2],s1,ls[dn],ms[dn],ls[n1],ms[n1],s2)
    # for n = 1:n_max for n2 = 1:n_max for dn = 1:n_max for n1 = 1:n_max ])
    # #
    # @time [
    #     M_component2(keff,Ns,n,dn,n1,n2,s1,s2)
    # for n = 1:n_max for dn = 1:n_max for n1 = 1:n_max for n2 = 1:n_max ];
    # #
    # @time [
    #      M_component(keff,Ns,ls[n],ms[n],ls[n2],ms[n2],s1,ls[dn],ms[dn],ls[n1],ms[n1],s2)
    # for n = 1:n_max for n2 = 1:n_max for dn = 1:n_max for n1 = 1:n_max ];
    return MM2
end
