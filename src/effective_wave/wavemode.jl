"calculate effective transmission angle θ_eff. We restrict -pi/2 < Re θ_eff < pi/2"
function transmission_angle(k::Complex{T}, k_eff::Complex{T}, θin; tol = 1e-8) where T<:Number
    snell(θ::Array{T}) = abs(k*sin(θin) - k_eff*sin(θ[1] + im*θ[2]))
    result = optimize(snell, [θin,0.]; x_tol= tol, g_tol= tol^2.0)

    if -pi/T(2) <= result.minimizer[1] <= pi/T(2)
        θ_eff = result.minimizer[1] + im*result.minimizer[2]
    else
        θ_eff = pi - result.minimizer[1] - im*result.minimizer[2]
    end
    return θ_eff
end

"The average effective transmitted scattering amplitudes for a given effective wavenumber k_eff. Assumes there exists only one k_eff.
The function returns an array A, where
AA(x,y,m,s) = im^m*exp(-im*m*θ_eff)*A[m + max_basis_order +1,s]*exp(im*k_eff*(cos(θ_eff)*x + sin(θin)*y))
where (x,y) are coordinates in the halfspace, m-th hankel order, s-th species,  and AA is the ensemble average scattering coefficient."
function effective_wavemodes(ω::T, k_eff::Complex{T}, medium::Acoustic{T,2}, species::Species{T};
        tol::T = 1e-5,
        kws...) where T<:Number

    MM = effectivewave_system(ω, medium, species; tol=tol, kws...)

    # calculate effective amplitudes
    MM_svd = svd(MM(k_eff))
    if MM_svd.S[end] > tol*T(100) # no analytic guarantee on this tolerance.
        @warn("The effective wavenumber gave an eigenvalue $(MM_svd.S[end]), which is larger than the tol = $( T(100)*tol). Note the dispersion equation is determined by the chosen medium and species.")
    end

    S = length(species)

    A_null = MM_svd.V[:,end] # eignvector of smallest eigenvalue
    A_null = reshape(A_null, (:,S)) # A_null[:,j] = [A[-ho,j],A[-ho+1,j],...,A[ho,j]]

    return A_null
end
