
function wavenumber(ωs::AbstractVector{T}, medium::Medium{T}, species::Vector{Specie{T}};
        wavenumber_initial_guess = :auto, kws...) where T<:Number
    ks = Vector{Complex{T}}(length(ωs))
    ks[1] = wavenumber(ωs[1], medium, species; wavenumber_initial_guess = wavenumber_initial_guess, kws...)
    ks[2] = wavenumber(ωs[2], medium, species; wavenumber_initial_guess = wavenumber_initial_guess, kws...)
    map(3:length(ks)) do i
        k0 = ks[i-2] + (ωs[i] - ωs[i-2])*(ks[i-1] - ks[i-2])/(ωs[i-1] - ωs[i-2])
        ks[i] = wavenumber(ωs[i], medium, species; wavenumber_initial_guess = k0, kws...)
    end
    return ks
end

function wavenumber(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; tol = 1e-8,
        wavenumber_initial_guess = :auto,
        hankel_order = maximum_hankel_order(ω, medium, species; tol=100*tol),
        radius_multiplier = 1.005,
        kws...) where T<:Number

    warn("This function is now depricated as it's too unstable. Use wavenumbers(ω::T, medium::Medium{T}, Vector{Specie{T}}), which will return many effective wavenumbers.")
    k = ω/medium.c
    S = length(species)
    ho = hankel_order
    Z_l_n = Zn_matrix(ω, medium, species; hankel_order = ho)

    r = maximum(s.r for s in species)
    φ = sum(volume_fraction.(species))
    if wavenumber_initial_guess == :auto
        if abs(r*k) > 0.2 && φ < 0.4
            # use low volume fraction (valid for any ω) effective wavenumber
            k0 = wavenumber_low_volfrac(ω, medium, species)
        else
            # use low frequency (valid for any φ) effective wavenumber
            eff_medium = effective_medium(medium, species)
            k0 = ω/eff_medium.c
        end
    else
        k0 = wavenumber_initial_guess
    end
    initial_k0 = [real(k0), 0.95*imag(k0)]

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    function M(keff,j,l,m,n)
        (n==m ? 1.0:0.0)*(j==l ? 1.0:0.0) + 2.0pi*species[l].num_density*Z_l_n[l,n]*
            Nn(n-m,k*as[j,l],keff*as[j,l])/(k^2.0-keff^2.0)
    end

    # this matrix is needed to calculate the eigenvectors
    MM(keff::Complex{T}) = reshape(
        [M(keff,j,l,m,n) for m in -ho:ho, j = 1:S, n in -ho:ho, l = 1:S]
    , ((2ho+1)*S, (2ho+1)*S))

    constraint(keff_vec::Array{T}) = ( (keff_vec[2] < zero(T)) ? one(T):zero(T))*(-1 + exp(-T(100.0)*keff_vec[2]))
    detMM2(keff_vec::Array{T}) =  constraint(keff_vec) + map(x -> real(x*conj(x)), det(MM(keff_vec[1]+im*keff_vec[2])))

    #Note that there is not a unique effective wavenumber. The root closest to k_eff = 0.0 + 0.0im seems to be the right one, the others lead to strange transmission angles and large amplitudes As.
    result = optimize(detMM2, initial_k0;  g_tol = tol^2.0, f_tol = tol^4.0)
    if !(result.f_converged || result.x_converged || result.g_converged) # if nothing converged
        result = optimize(detMM2, [0.,0.];  g_tol = tol^2.0, f_tol = tol^4.0)
        if !(result.f_converged || result.x_converged || result.g_converged)
            warn("Local optimisation did not converge")
        end
    end
    # Check result
    k_eff = result.minimizer[1] + im*result.minimizer[2]
    # in case wave travelling opposite direction was found
    if imag(k_eff) < zero(T) k_eff = - k_eff end
    MM_svd = svd(MM(k_eff))
    if last(MM_svd[2]) > T(4)*tol
        warn("Local optimisation was unsucessful at finding an effective wavenumber for ω = $ω and max(radius)*k = $(abs(r*k)).")
        warn("$(last(MM_svd[2])) was the smallest eigenvalue value (should be zero) of the effective wavenumber matrix equation.")
    end

    return k_eff
end
