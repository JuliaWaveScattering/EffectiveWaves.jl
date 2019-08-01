function dispersion_function(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        dim = 2, tol::T = 1e-4,
        kws...) where T<:Number

    low_tol = max(1e-4, tol) # a tolerance used for a first pass with time_limit

    MM = if dim == 2
        wavematrix2D(ω, medium, species; tol=tol, kws...)
    elseif dim == 3
        wavematrix3D(ω, medium, species; tol=tol, kws...)
    elseif dim == :plane
        wavematrix3DPlane(ω, medium, species; tol=tol, kws...)
    else error("the dimension dim = $dim has no dispersion function implemented.")
    end

    # the constraint uses keff_vec[2] < -low_tol to better specify solutions where imag(k_effs)<0
    constraint(keff_vec::Array{T}) = ( (keff_vec[2] < -low_tol) ? one(T) : zero(T))*(-1 + exp(-T(100.0)*keff_vec[2]))
    detMM2(keff_vec::Array{T}) =  constraint(keff_vec) +
        map(x -> real(x*conj(x)), det(MM(keff_vec[1]+im*keff_vec[2])))

    return detMM2
end

function wavematrix2D(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        dim = 2, tol::T = 1e-4,
        hankel_order::Int = maximum_hankel_order(ω, medium, species; tol=tol),
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; hankel_order = hankel_order, dim = dim),
        kws...) where T<:Number

    k = real(ω/medium.c)
    S = length(species)
    ho = hankel_order

    len = (2ho+1) * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    function M_component(keff,j,l,m,n)
        (n == m ? 1.0 : 0.0)*(j == l ? 1.0 : 0.0) - 2.0pi*species[l].num_density*t_vecs[l][m+ho+1]*
            kernelN(n-m,k*as[j,l],keff*as[j,l])/(k^2.0-keff^2.0)
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

function wavematrix3D(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        dim = 3, tol::T = 1e-4,
        hankel_order::Int = maximum_hankel_order(ω, medium, species; tol=tol),
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; hankel_order = hankel_order, dim = dim),
        kws...) where T<:Number

    k = real(ω/medium.c)
    S = length(species)
    ho = hankel_order
    len = (ho+1)^4 * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    function M_component(keff,l,m,l2,m2,s1,dl,dm,l1,m1,s2)::Complex{T}
        (m == dm && l == dl && m1 == m2 && l1 == l2 && s1 == s2 ? 1.0 : 0.0) +
        as[s1,s2] * species[s2].num_density * t_vecs[s1][l+ho+1] *
        sum(l3 ->
            if abs(m1-m2) <= l3
                gs = gaunt_coefficients(l,m,dl,dm,l3,m1-m2) *
                gaunt_coefficients(l1,m1,l2,m2,l3,m1-m2)
                if abs(gs) < eps(T)
                    gs * kernelN(l3,k*as[s1,s2],keff*as[s1,s2]; dim = dim)
                else zero(Complex{T})
                end
            else
                zero(Complex{T})
            end
        , 0:ho) / (keff^2.0 - k^2.0)
    end

    # The order of the indices below has not been carefully planned
    function MM(keff::Complex{T})::Matrix{Complex{T}}
        ind2 = 1
        for s2 = 1:S for dl = 0:ho for dm = -dl:dl for l1 = 0:ho for m1 = -l1:l1
            ind1 = 1
            for s1 = 1:S for l = 0:ho for m = -l:l for l2 = 0:ho for m2 = -l2:l2
                MM_mat[ind1, ind2] = M_component(keff,l,m,l2,m2,s1,dl,dm,l1,m1,s2)
                ind1 += 1
            end end end end end
            ind2 += 1
        end end end end end
        return MM_mat
    end

    return MM
end

function wavematrix3D_allocate(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        dim = 3, tol::T = 1e-4,
        hankel_order::Int = maximum_hankel_order(ω, medium, species; tol=tol),
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; hankel_order = hankel_order, dim = dim),
        kws...) where T<:Number

    k = real(ω/medium.c)
    S = length(species)
    ho = hankel_order
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

    function M_component2(keff,n,dn,n1,n2,s1,s2)::Complex{T}
        (n == dn && n1 == n2 && s1 == s2 ? 1.0 : 0.0) +
        as[s1,s2] * species[s2].num_density * t_vecs[s1][ls[n]+ho+1] *
        sum(n3 ->
        begin
            cm = cs[n3,n,dn] * cs[n3,n1,n2]
            if abs(cm) < eps(T)
                zero(Complex{T})
            else
                cm * kernelN(ls[n3],k*as[s1,s2],keff*as[s1,s2]; dim = dim)
            end
        end
        , 1:n_max) / (keff^2.0 - k^2.0)
    end
    function MM2(keff::Complex{T})
        ind2 = 1
        for s2 = 1:S for dn = 1:n_max for n1 = 1:n_max
            ind1 = 1
            for s1 = 1:S for n = 1:n_max for n2 = 1:n_max
                MM_mat[ind1, ind2] = M_component2(keff,n,dn,n1,n2,s1,s2)
                ind1 += 1
            end end end
            ind2 += 1
        end end end
        return MM_mat
    end
    # keff = 1.0 + rand()*im
    # s1 = s2 = 1

    # norm([
    #     M_component2(keff,n,dn,n1,n2,s1,s2) - M_component(keff,ls[n],ms[n],ls[n2],ms[n2],s1,ls[dn],ms[dn],ls[n1],ms[n1],s2)
    # for n = 1:n_max for n2 = 1:n_max for dn = 1:n_max for n1 = 1:n_max ])
    # #
    # @time [
    #     M_component(keff,n,dn,n1,n2,s1,s2)
    # for n = 1:n_max for dn = 1:n_max for n1 = 1:n_max for n2 = 1:n_max ];
    #
    # @time [
    #      M_component(keff,ls[n],ms[n],ls[n2],ms[n2],s1,ls[dn],ms[dn],ls[n1],ms[n1],s2)
    # for n = 1:n_max for n2 = 1:n_max for dn = 1:n_max for n1 = 1:n_max ];
    return MM2
end


function wavematrix3DPlane(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        θ_inc::T = T(0), # incident on the halfspace z >0, with kp_vec being in the y-z plane
        φ_inc::T = T(0),
        dim = 3, tol::T = 1e-4,
        hankel_order::Int = maximum_hankel_order(ω, medium, species; tol=tol),
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; hankel_order = hankel_order, dim = dim),
        kws...) where T<:Number

    k = real(ω/medium.c)
    S = length(species)
    ho = hankel_order
    len = (ho+1)^2 * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    Ys = spherical_harmonics(ho, θ_inc, φ_inc);
    ls, ms = spherical_harmonics_indices(ho);

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    function M_component(keff,l,m,s1,dl,dm,s2)
        (m == dm && l == dl && s1 == s2 ? 1.0 : 0.0) +
        T(4pi) * as[s1,s2] * species[s2].num_density * t_vecs[s1][l+ho+1] *
        sum(
            if ms[n1] != dm - m
                zero(T)
            else
                Complex{T}(im)^(-ls[n1]) * Ys[n1] * gaunt_coefficients(dl,dm,l,m,ls[n1],ms[n1]) *
                kernelN(ls[n1],k*as[s1,s2],keff*as[s1,s2]; dim = dim)
            end
        for n1 in eachindex(Ys)) / (keff^2.0 - k^2.0)
    end

    function MM(keff::Complex{T})
        ind2 = 1
        for s2 = 1:S, dl = 0:ho for dm = -dl:dl
            ind1 = 1
            for s1 = 1:S, l = 0:ho for m = -l:l
                MM_mat[ind1, ind2] = M_component(keff,l,m,s1,dl,dm,s2)
                ind1 += 1
            end end
            ind2 += 1
        end end
        return MM_mat
    end

    return MM
end
