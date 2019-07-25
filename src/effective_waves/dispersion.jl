function dispersion_function(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        dim = 2, tol::T = 1e-4,
        kws...) where T<:Number

    low_tol = max(1e-4, tol) # a tolerance used for a first pass with time_limit

    MM = if dim == 2
        wavematrix2D(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; tol=tol, kws...)
    else
        wavematrix3D(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; tol=tol, kws...)
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
            end
            end
            ind2 += 1
        end
        end
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
    function M_component(keff,l,m,l2,m2,s2,dl,dm,l1,m1,s1)
        (m == dm && l == dl ? 1.0 : 0.0)*(m1 == m2 && l1 == l2 ? 1.0 : 0.0) +
<<<<<<< HEAD
        as[s2,s1] * species[s1].num_density * t_vecs[s2][l+ho+1] *
        sum(l3 ->
            if abs(m1-m2) <= l3
                gaunt_coefficients(l, m, dl,dm,l3,m1-m2) *
                gaunt_coefficients(l1,m1,l2,m2,l3,m1-m2) *
                kernelN(l3,k*as[s2,s1],keff*as[s2,s1]; dim = dim)
            else
                zero(T)
            end
        , 0:ho) / (keff^2.0 - k^2.0)
    end

    # The order of the indices below has not been carefully planned
    function MM(keff::Complex{T})
        ind2 = 1
        for s1 = 1:S, dl = 0:ho, l1 = 0:ho for dm = -dl:dl, m1 = -l1:l1
            ind1 = 1
            for s2 = 1:S, l = 0:ho, l2 = 0:ho for m = -l:l, m2 = -l2:l2
                MM_mat[ind1, ind2] = M_component(keff,l,m,l2,m2,s2,dl,dm,l1,m1,s1)
                ind1 += 1
            end end
            ind2 += 1
        end end
        return MM_mat
    end

    return MM
end


function wavematrix3DPlane(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
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
    function M_component(keff,l,m,l2,m2,s2,dl,dm,l1,m1,s1)
        (m == dm && l == dl ? 1.0 : 0.0)*(m1 == m2 && l1 == l2 ? 1.0 : 0.0) +
=======
>>>>>>> 2820b780161c36e6c357cdb94c56194e607ad592
        as[s2,s1] * species[s1].num_density * t_vecs[s2][l+ho+1] *
        sum(l3 -> sum(m3 ->
            gaunt_coefficients(l, m, dl,dm,l3,m3) *
            gaunt_coefficients(l1,m1,l2,m2,l3,m3) *
            kernelN(l3,k*as[s2,s1],keff*as[s2,s1]; dim = dim)
        , -l3:l3), 0:ho) / (keff^2.0 - k^2.0)
    end

    # The order of the indices below has not been carefully planned
    function MM(keff::Complex{T})
        ind2 = 1
        for s1 = 1:S, dl = 0:ho, l1 = 0:ho for dm = -dl:dl, m1 = -l1:l1
            ind1 = 1
            for s2 = 1:S, l = 0:ho, l2 = 0:ho for m = -l:l, m2 = -l2:l2
                MM_mat[ind1, ind2] = M_component(keff,l,m,l2,m2,s2,dl,dm,l1,m1,s1)
                ind1 += 1
            end end
            ind2 += 1
        end end
        return MM_mat
    end

    return MM
end
