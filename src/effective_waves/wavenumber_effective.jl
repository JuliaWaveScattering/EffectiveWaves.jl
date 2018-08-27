# Here we calculate the effective wavenumber and effective wave amplitudes, without any restriction to the volume fraction of particles and incident wave frequency.

Nn(n,x,y) = x*diffhankelh1(n,x)*besselj(n,y) - y*hankelh1(n,x)*diffbesselj(n,y)

# include depricated function to find a single effective wavenumber, when in fact there are many. The code is still used in tests and gives many correct results
include("wavenumber_single.jl")

" Returns all the transmitted effective wavenumbers"
wavenumbers(ω::T, medium::Medium{T}, specie::Specie{T}; kws...) where T<:Number = wavenumbers(ω, medium, [specie]; kws...)

function wavenumbers(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; tol::T = 1e-6,
        hankel_order::Int = maximum_hankel_order(ω, medium, species; tol=tol),
        mesh_points::Int = 5, mesh_size::T = 0.5,
        max_Imk::T = 0.0, max_Rek::T = 0.0,
        time_limit::T = 1.0,
        radius_multiplier::T = 1.005,
        kws...) where T<:Number

    k = ω/medium.c
    S = length(species)
    ho = hankel_order

    Z_l_n = Zn_matrix(ω, medium, species; hankel_order = ho)

    r = maximum(s.r for s in species)
    φ = sum(volume_fraction.(species))

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    function M_component(keff,j,l,m,n)
        (n==m ? 1.0:0.0)*(j==l ? 1.0:0.0) + 2.0pi*species[l].num_density*Z_l_n[l,n]*
            Nn(n-m,k*as[j,l],keff*as[j,l])/(k^2.0-keff^2.0)
    end

    # this matrix is needed to calculate the eigenvectors
    MM(keff::Complex{T}) = reshape(
        [M_component(keff,j,l,m,n) for m in -ho:ho, j = 1:S, n in -ho:ho, l = 1:S]
    , ((2ho+1)*S, (2ho+1)*S))

    constraint(keff_vec::Array{T}) = ( (keff_vec[2] < zero(T)) ? one(T):zero(T))*(-1 + exp(-T(100.0)*keff_vec[2]))
    detMM2(keff_vec::Array{T}) =  constraint(keff_vec) + map(x -> real(x*conj(x)), det(MM(keff_vec[1]+im*keff_vec[2])))

    kφ = wavenumber_low_volfrac(ω, medium, species; verbose = false)
    eff_medium = effective_medium(medium, species)
    k0 = ω/eff_medium.c
    if isnan(k0) k0 = kφ end

    if max_Rek == 0.0
        dk_x = max(real(k0),abs(real(kφ))) * mesh_size
        max_Rek = mesh_points * dk_x
    else dk_x = max_Rek/mesh_points
    end
    if max_Imk == 0.0
        dk_y = abs(imag(kφ)) * mesh_size
        max_Imk = mesh_points * dk_y
    else dk_y = max_Imk/mesh_points
    end

    kx = -max_Rek:dk_x:max_Rek
    ky = 0.0:dk_y:max_Imk
    kins = [[x,y] for x in kx, y in ky]
    #Note that there is not a unique effective wavenumber. The root closest to k_eff = 0.0 + 0.0im seems to be the right one, the others lead to strange transmission angles and large amplitudes As.
    k_vecs = map(kins) do kin
        result = optimize(detMM2, kin; time_limit = time_limit)
        if result.minimum < 100*tol
            result.minimizer
        else
            [0.,-1.0]
        end
    end
    k_vecs = sort(k_vecs[:]; by = x -> x[2])

    digs = -3 - Int(round(log(10,tol))) # number of digits to round to check if equal
    k_vecs = map(groupby(k_vec -> round.(k_vec,digs), k_vecs)) do g
        res = mean(g)
    end
    deleteat!(k_vecs, find(k_vec[2] < 0 for k_vec in k_vecs))

    # Here we refine the effective wavenumbers
    k_vecs = map(k_vecs) do k_vec
        res = optimize(detMM2, k_vec;  g_tol = tol^2.0, f_tol = tol^4.0)
        res.minimizer
    end

    k_vecs = sort(k_vecs[:]; by = x -> x[2])
    k_effs = map(groupby(k_vec -> round.(k_vec,digs), k_vecs)) do g
        res = mean(g)
        res[1] + res[2]*im
    end

    return k_effs
end

# detMM2(keff_vec::Array{T}) =  map(x -> real(x*conj(x)), det(MM(keff_vec[1]+im*keff_vec[2])))

# function detMM!(F,x)
#     F[1] = abs(det(MM(x[1]+im*x[2])))
# end

# Alternative solvers
# res = nlsolve(detMM!,initial_k_eff; iterations = 10000, factor=2.0)
# k_eff_nl = res.zero[1] + im*res.zero[2]
# lower = [0.,0.]; upper = [T(2)*k0,k0]
# result = optimize(detMM2, initial_k0, lower, upper; g_tol = tol^2.0, f_tol = tol^4.0)
