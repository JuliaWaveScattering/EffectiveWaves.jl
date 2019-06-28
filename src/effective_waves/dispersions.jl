# Here are the available dispersion equations


function dispersion_2d_plane(ω::T, species::Vector{Specie{T}}, medium::Medium{T}, species::Vector{Specie{T}};
        tol::T = 1e-6,
        hankel_order::Int = maximum_hankel_order(ω, medium, species; tol=tol),
        radius_multiplier::T = T(1.005),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, species; hankel_order = hankel_order)
    ) where T<:Number

    # low_tol = max(1e-4, tol) # a tolerance used for a first pass with time_limit
    # k = real(ω/medium.c)
    # S = length(species)
    # ho = hankel_order
    #
    # as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    # function M_component(keff,j,l,m,n)
    #     (n == m ? 1.0 : 0.0)*(j == l ? 1.0 : 0.0) - 2.0pi*species[l].num_density*t_vecs[l][m+ho+1]*
    #         Nn(n-m,k*as[j,l],keff*as[j,l])/(k^2.0-keff^2.0)
    # end
    #
    # # this matrix is needed to calculate the eigenvectors
    # MM(keff::Complex{T}) = reshape(
    #     [M_component(keff,j,l,m,n) for m in -ho:ho, j = 1:S, n in -ho:ho, l = 1:S]
    # , ((2ho+1)*S, (2ho+1)*S))
    #
    # # the constraint uses keff_vec[2] < -low_tol to better specify solutions where imag(k_effs)<0
    # constraint(keff_vec::Array{T}) = ( (keff_vec[2] < -low_tol) ? one(T) : zero(T))*(-1 + exp(-T(100.0)*keff_vec[2]))
    # detMM2(keff_vec::Array{T}) =  constraint(keff_vec) + map(x -> real(x*conj(x)), det(MM(keff_vec[1]+im*keff_vec[2])))
end
