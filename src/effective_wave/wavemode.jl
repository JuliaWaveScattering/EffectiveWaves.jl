
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

"returns a number a, such that a*As_eff will cancel an incident wave plane wave with incident angle θin."
function scale_amplitudes_effective(ω::T, wave_eff::EffectiveWave{T},
        medium::Acoustic{T,2}, species::Species{T};
        θin::T = 0.0, tol = 1e-8) where T<:Number

    k = ω/medium.c
    θ_eff = wave_eff.θ_eff
    ho = wave_eff.basis_order
    amps = wave_eff.amplitudes
    S = length(species)

    sumAs = T(2)*sum(
            exp(im*n*(θin - θ_eff))*number_density(species[l])*amps[n+ho+1,l]
    for n = -ho:ho, l = 1:S)
    a = im*k*cos(θin)*(wave_eff.k_eff*cos(θ_eff) - k*cos(θin))/sumAs

    return a
end

function wienerhopf_wavemodes(ω::T, k_eff::Complex{T}, medium::Acoustic{T,2}, species::Species{T}; kws...) where T<:AbstractFloat
    return wienerhopf_wavemodes(ω, [k_eff], medium, species; kws...)
end


"The average effective transmitted wavemodes according to the Wiener-Hopf method.
The function returns an array A, where
AA(x,y,0,1) = A[1,1]*exp(im*k_eff*(cos(θ_eff)*x + sin(θin)*y))
where (x,y) are coordinates in the halfspace  and AA is the ensemble average scattering coefficient. Method currently only implemented for 1 species and for monopole scatterers."
function wienerhopf_wavemodes(ω::T, k_effs::Vector{Complex{T}}, medium::Acoustic{T,2}, species::Species{T};
        tol::T = 1e-6, θin::T = 0.0,
        basis_order::Int = 0,
        num_coefs::Int = 10000,
        maxZ::T = T(100)*maximum(outer_radius(s) * s.exclusion_distance for s in species) + T(100),
        kws...
    ) where T<:AbstractFloat

    k = ω/medium.c
    ho = basis_order

    if ho > 0
        error("the Wiener Hopf method has not been implemented for `basis_order` = $(basis_order). Method currently only works for basis_order = 0, i.e. monopole scatterers. ")
    end

    t_vecs = get_t_matrices(medium, species, ω, basis_order)

    # t_vecs = t_vectors(ω, medium, species; basis_order = ho)
    as = [
        (outer_radius(s1) * s1.exclusion_distance + outer_radius(s2) * s2.exclusion_distance)
    for s1 in species, s2 in species]

    # differentiate in S ( = a[j,l] k_eff ) and evaluate at S
    function dSQ0_eff(S,j,l)
        m = 0

        T(2) / (S^T(2) - (k*as[j,l])^T(2)) * (
            S +
            pi * as[j,l]^T(2) * number_density(species[l]) * t_vecs[l][m+ho+1,m+ho+1] * (
                k*as[j,l]*diffhankelh1(0,k*as[j,l])*diffbesselj(0,S) -
                hankelh1(0,k*as[j,l])*diffbesselj(0,S) -
                S*hankelh1(0,k*as[j,l])*diffbesselj(0,S,2)
            )
        )
    end

    # kernelN(0,k*a12,Z) = k*a12*diffhankelh1(0,k*a12)*besselj(0,Z) - Z*hankelh1(0,k*a12)*diffbesselj(0,Z)
    # DZkernelN(0,k*a12,Z) = k*a12*diffhankelh1(0,k*a12)*diffbesselj(0,Z) - Z*hankelh1(0,k*a12)*diffbesselj(0,Z,2) -

    sToS(s,j::Int,l::Int) = (real(s) >= 0) ? sqrt(s^2 + (k*as[j,l]*sin(θin))^2) : -sqrt(s^2 + (k*as[j,l]*sin(θin))^2)

    function q(s,j,l,m,n)
        (n == m ? T(1) : T(0)) * (j == l ? T(1) : T(0)) +
        T(2) * pi * as[j,l]^T(2) * number_density(species[l]) * t_vecs[l][m+ho+1,m+ho+1] *
        kernelN(n-m, k*as[j,l], sToS(s,j,l)) / (s^T(2) - (k*as[j,l]*cos(θin))^T(2))
    end

    Zs = LinRange(T(100),1/(10*tol),3000)
    maxZ = Zs[findfirst(Z -> abs(log(q(Z,1,1,0,0))) < 10*tol, Zs)]

    function Ψp(s, maxZ::T = maxZ, num_coefs::Int = num_coefs)
        Q(z) = log(q(z,1,1,0,0))/(z - s)
        xp = as[1,1]*k*cos(θin)*(-1.0+1.0im)
        (s + k*as[1,1]*cos(θin)) * exp(
            (T(1.0)/(T(2)*pi*im)) * (
                sum(Fun(Q, Segment(-maxZ,xp), num_coefs)) +
                sum(Fun(Q, Segment(xp,-xp), num_coefs)) +
                sum(Fun(Q, Segment(-xp,maxZ), num_coefs))
            )
        )
    end

    Ψp_a = Ψp(k*as[1,1]*cos(θin))
    # has been tested against Mathematica for at least one k_eff
    dSΨ00(S) = (S^2 - (k*as[1,1])^2)*dSQ0_eff(S,1,1) # + T(2) * S * Q0(S,1,1,0,0)
    # last term left out becuase Q0(k_eff,1,1,0,0) = 0.

    return map(k_effs) do k_eff
        θ_eff = transmission_angle(k, k_eff, θin; tol = tol)
        Ψp_eff = Ψp(k_eff*as[1,1]*cos(θ_eff))
        t_vecs[1][ho+1,ho+1] * T(2)*k*as[1,1]*(cos(θin)/cos(θ_eff))*(Ψp_eff/Ψp_a)/dSΨ00(k_eff*as[1,1])
    end
end
