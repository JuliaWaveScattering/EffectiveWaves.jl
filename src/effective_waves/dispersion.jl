function dispersion_function(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        dim = 2, tol::T = 1e-4, symmetry = :plane,
        kws...) where T<:Number

    low_tol = max(1e-4, tol) # a tolerance used for a first pass with time_limit

    MM = if dim == 2 && symmetry == :plane
        wavematrix2D(ω, medium, species; tol=tol, kws...)
    elseif dim == 3 && symmetry == :plane
        wavematrix3DPlane(ω, medium, species; tol=tol, kws...)
    elseif dim == 3 && symmetry == :azimuth
        wavematrix3DAzimuth(ω, medium, species; tol=tol, kws...)
    elseif dim == 3 && symmetry == :planeazimuth
        wavematrix3DPlaneAzimuth(ω, medium, species; tol=tol, kws...)
    elseif dim == 3 && symmetry == :none
        wavematrix3D(ω, medium, species; tol=tol, kws...)
    else error("the dimension dim = $dim has no dispersion function implemented.")
    end

    # the constraint uses keff_vec[2] < -low_tol to better specify solutions where imag(k_effs)<0
    constraint(keff_vec::Array{T}) = (keff_vec[2] < -low_tol) ? (-one(T) + exp(-T(100.0)*keff_vec[2])) : zero(T)

    detMM(keff) = det(MM(keff))
    detMM2 = if dim == 3 && symmetry != :plane
        function(keff_vec::Array{T})
            constraint(keff_vec) + sqrt(abs(detMM(keff_vec[1]+im*keff_vec[2])))
        end
    else
        function detMM2(keff_vec::Array{T})
            constraint(keff_vec) + abs2(detMM(keff_vec[1]+im*keff_vec[2]))
        end
    end

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
        hankel_order::Int = 2,
        hankel_order_field::Int = 2*hankel_order,
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; hankel_order = hankel_order, dim = dim),
        kws...) where T<:Number

    k = real(ω/medium.c)
    S = length(species)
    L = hankel_order
    L1 = hankel_order_field
    len = (L1+1)^2 * (L+1)^2 * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    function M_component(keff,Ns,l,m,l2,m2,s1,dl,dm,l1,m1,s2)::Complex{T}
        (m == dm && l == dl && m1 == m2 && l1 == l2 && s1 == s2 ? 1.0 : 0.0) +
        as[s1,s2] * species[s2].num_density * t_vecs[s1][l+L+1] *
        sum(l3 ->
            if abs(m1-m2) <= l3
                gaunt_coefficients(l,m,dl,dm,l3,m1-m2) *
                gaunt_coefficients(l1,m1,l2,m2,l3,m1-m2) * Ns[l3+1,s1,s2]
            else
                zero(Complex{T})
            end
        , 0:L1) / (keff^2.0 - k^2.0)
    end

    # The order of the indices below is important
    function MM(keff::Complex{T})::Matrix{Complex{T}}
        Ns = [kernelN(l3,k*as[s1,s2],keff*as[s1,s2]; dim = dim) for l3 = 0:L1, s1 = 1:S, s2 = 1:S]
        ind2 = 1
        for s2 = 1:S for dl = 0:L for dm = -dl:dl for l1 = 0:L1 for m1 = -l1:l1
            ind1 = 1
            for s1 = 1:S for l = 0:L for m = -l:l for l2 = 0:L1 for m2 = -l2:l2
                MM_mat[ind1, ind2] = M_component(keff,Ns,l,m,l2,m2,s1,dl,dm,l1,m1,s2)
                ind1 += 1
            end end end end end
            ind2 += 1
        end end end end end
        return MM_mat
    end

    return MM
end

function wavematrix3DAzimuth(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        dim = 3, tol::T = 1e-4,
        hankel_order::Int = 1,
        hankel_order_field::Int = 2*hankel_order,
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; hankel_order = hankel_order, dim = dim),
        kws...) where T<:Number

    k = real(ω/medium.c)
    S = length(species)
    L = hankel_order
    L1 = hankel_order_field

    len = Int(1 - L*(2 + L)*(L - 3*L1 - 2)/3 + L1)

    MM_mat = Matrix{Complex{T}}(undef,len,len)

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    # the index for the T-matrix below needs to be changed when seperating correctly the 2D and 3D case.
    function M_component(keff,Ns,l,m,l2,s1,dl,dm,l1,s2)::Complex{T}
        minl3 = max(abs(m-dm),abs(dl-l),abs(l1-l2))
        maxl3 = min(dl+l,l1+l2)

        (m == dm && l == dl && l1 == l2 && s1 == s2 ? 1.0 : 0.0) +
        if abs(dm) <= min(l1) && abs(m) <= min(l2) && minl3 <= maxl3
            as[s1,s2] * species[s2].num_density * t_vecs[s1][l+L+1] *
            sum(l3 ->
                gaunt_coefficients(l,m,dl,dm,l3,m-dm) *
                gaunt_coefficients(l1,-dm,l2,-m,l3,m-dm) * Ns[l3+1,s1,s2]
            , minl3:maxl3) / (keff^2.0 - k^2.0)
        else
            zero(Complex{T})
        end
    end

    # The order of the indices below is important
    function MM(keff::Complex{T})::Matrix{Complex{T}}
        Ns = [kernelN(l3,k*as[s1,s2],keff*as[s1,s2]; dim = dim) for l3 = 0:min(2L1,2L), s1 = 1:S, s2 = 1:S]
        ind2 = 1
        for s2 = 1:S for dl = 0:L for dm = -dl:dl for l1 = abs(dm):L1
            ind1 = 1
            for s1 = 1:S for l = 0:L for m = -l:l for l2 = abs(m):L1
                MM_mat[ind1, ind2] = M_component(keff,Ns,l,m,l2,s1,dl,dm,l1,s2)
                ind1 += 1
            end end end end
            ind2 += 1
        end end end end
        return MM_mat
    end

    return MM
end

function wavematrix3DPlane(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        θ_inc::T = T(0), # incident on the halfspace z >0, with kp_vec being in the y-z plane
        φ_inc::T = T(0),
        dim = 3, tol::T = 1e-5,
        hankel_order::Int = 1,
        hankel_order_field::Int = 2*hankel_order,
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; hankel_order = hankel_order, dim = dim),
        kws...) where T<:Number

    k = real(ω/medium.c)
    S = length(species)
    ho = hankel_order
    L1 = hankel_order_field
    len = (ho+1)^2 * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    Ys = spherical_harmonics(L1, θ_inc, φ_inc);
    lm_to_n = lm_to_spherical_harmonic_index

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    function M_component(keff::Complex{T},Ns::Array{Complex{T}},l,m,s1,dl,dm,s2)::Complex{T}
        (m == dm && l == dl && s1 == s2 ? one(Complex{T}) : zero(Complex{T})) +
        Complex{T}(4pi) * as[s1,s2] * species[s2].num_density * t_vecs[s1][l+ho+1] *
        sum(
            Complex{T}(im)^(-l1) * Ys[lm_to_n(l1,dm-m)] * Ns[l1+1,s1,s2] *
            gaunt_coefficients(dl,dm,l,m,l1,dm-m)
        for l1 in max(abs(dm-m),abs(dl-l)):(dl+l)) / (keff^2.0 - k^2.0)
    end

    function MM(keff::Complex{T})::Matrix{Complex{T}}
        Ns = [kernelN(l1,k*as[s1,s2],keff*as[s1,s2]; dim = dim) for l1 = 0:L1, s1 = 1:S, s2 = 1:S]
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

function wavematrix3DPlaneAzimuth(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        dim = 3, tol::T = 1e-4,
        hankel_order::Int = 1,
        hankel_order_field::Int = 2*hankel_order,
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; hankel_order = hankel_order, dim = dim),
        kws...) where T<:Number

    k = real(ω/medium.c)
    S = length(species)
    ho = hankel_order
    L1 = hankel_order_field
    len = (ho+1) * S
    MM_mat = Matrix{Complex{T}}(undef,len,len)

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    function M_component(keff::Complex{T},Ns::Array{Complex{T}},l,s1,dl,s2)::Complex{T}
        (l == dl && s1 == s2 ? one(Complex{T}) : zero(Complex{T})) +
        Complex{T}(4pi) * as[s1,s2] * species[s2].num_density * t_vecs[s1][l+ho+1] *
        sum(
            Complex{T}(im)^(-l1) * sqrt((2*l1+1)/(4pi) ) * Ns[l1+1,s1,s2] *
            gaunt_coefficients(dl,0,l,0,l1,0)
        for l1 in abs(dl-l):(dl+l)) / (keff^2.0 - k^2.0)
    end

    function MM(keff::Complex{T})::Matrix{Complex{T}}
        Ns = [kernelN(l1,k*as[s1,s2],keff*as[s1,s2]; dim = dim) for l1 = 0:L1, s1 = 1:S, s2 = 1:S]
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
