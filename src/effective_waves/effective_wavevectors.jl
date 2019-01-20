
"The average effective transmitted scattering amplitudes for a given effective wavenumber k_eff. Assumes there exists only one k_eff.
The function returns an array A, where
AA(x,y,m,s) = im^m*exp(-im*m*θ_eff)*A[m + max_hankel_order +1,s]*exp(im*k_eff*(cos(θ_eff)*x + sin(θin)*y))
where (x,y) are coordinates in the halfspace, m-th hankel order, s-th species,  and AA is the ensemble average scattering coefficient."
function effective_wavevectors(ω::T, k_eff::Complex{T}, medium::Medium{T}, species::Vector{Specie{T}};
        tol::T = 1e-5,
        radius_multiplier::T = 1.005,
        hankel_order::Int = maximum_hankel_order(ω, medium, species; tol=tol)
        ) where T<:Number

    k = ω/medium.c
    S = length(species)
    ho = hankel_order

    t_vecs = t_vectors(ω, medium, species; hankel_order = hankel_order)

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    function MM(keff,j,l,m,n)
        - (n == m ? 1.0 : 0.0)*(j == l ? 1.0 : 0.0) + 2.0pi*species[l].num_density*t_vecs[l][m+ho+1]*
            Nn(n-m,k*as[j,l],keff*as[j,l])/(k^2.0-keff^2.0)
    end

    # calculate effective amplitudes
    MM_svd = svd(
        reshape(
            [MM(k_eff,j,l,m,n) for m in -ho:ho, j in 1:S, n in -ho:ho, l in 1:S]
        , ((2ho+1)*S, (2ho+1)*S))
    )
    if MM_svd.S[end] > tol*T(40) # no analytic guarantee on this tolerance.
        @warn("The effective wavenumber gave an eigenvalue $(MM_svd.S[end]), which is larger than the tol = $( T(40)*tol). Note the secular equation is determined by the chosen medium and species.")
    end

    A_null = MM_svd.V[:,end] # eignvector of smallest eigenvalue
    A_null = reshape(A_null, (2*ho+1,S)) # A_null[:,j] = [A[-ho,j],A[-ho+1,j],...,A[ho,j]]

    return A_null
end

"returns a number a, such that a*As_eff will cancel an incident wave plane wave with incident angle θin."
function scale_amplitudes_effective(ω::T, wave_eff::EffectiveWave{T},
        medium::Medium{T}, species::Vector{Specie{T}};
        θin::T = 0.0, tol = 1e-8) where T<:Number

    k = ω/medium.c
    θ_eff = wave_eff.θ_eff
    ho = wave_eff.hankel_order
    amps = wave_eff.amplitudes
    S = length(species)

    sumAs = T(2)*sum(
            exp(im*n*(θin - θ_eff))*species[l].num_density*amps[n+ho+1,l]
    for n = -ho:ho, l = 1:S)
    a = im*k*cos(θin)*(wave_eff.k_eff*cos(θ_eff) - k*cos(θin))/sumAs

    return a
end

"The average effective transmitted wavevectors according to the Wiener-Hopf method.
The function returns an array A, where
AA(x,y,0,1) = A[1,1]*exp(im*k_eff*(cos(θ_eff)*x + sin(θin)*y))
where (x,y) are coordinates in the halfspace  and AA is the ensemble average scattering coefficient. Method currently only implemented for 1 species and for monopole scatterers."
function wienerhopf_wavevectors(ω::T, k_eff::Complex{T}, medium::Medium{T}, species::Vector{Specie{T}};
        tol::T = 1e-6, θin::T = 0.0,
        radius_multiplier::T = 1.005,
        hankel_order::Int = 0,
        num_coefs::Int = 10000
        ) where T<:AbstractFloat

    k = ω/medium.c
    ho = hankel_order
    θ_eff = transmission_angle(k, k_eff, θin; tol = tol)

    t_vecs = t_vectors(ω, medium, species; hankel_order = ho)
    # Nn(0,k*a12,Z) = k*a12*diffhankelh1(0,k*a12)*besselj(0,Z) - Z*hankelh1(0,k*a12)*diffbesselj(0,Z)
    # Q(S) = Q0(S/as[1,1],1,1,0,0) = 1 + 2.0pi*species[l].num_density*as[1,1]^2*t_vecs[l][m+ho+1] * Nn(0,k*species[1].r,S) / (S^2 - k^2 * as[1,1]^2);

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]

    function Q0(keff,j,l,m,n)
        (n == m ? T(1) : T(0))*(j == l ? T(1) : T(0)) + T(2) * pi*species[l].num_density*t_vecs[l][m+ho+1]*
            Nn(n-m,k*as[j,l],keff*as[j,l])/(keff^T(2) - k^T(2))
    end

    # differentiate in S ( = a[j,l] k_eff )
    function dQ0(j,l)
        m = 0
        T(2) * pi*species[l].num_density*t_vecs[l][m+ho+1]*(-T(2)*k_eff)/(k_eff^T(2) - k^T(2))^2 *
        Nn(0,k*as[j,l],k_eff*as[j,l])/as[j,l] +
        T(2) * pi*species[l].num_density*t_vecs[l][m+ho+1]/(k_eff^T(2) - k^T(2))*(
            k*as[j,l]*diffhankelh1(0,k*as[j,l])*diffbesselj(0,k_eff*as[j,l]) -
            hankelh1(0,k*as[j,l])*diffbesselj(0,k_eff*as[j,l]) -
            k_eff*as[j,l]*hankelh1(0,k*as[j,l])*diffbesselj(0,k_eff*as[j,l],2)
        )
    end

    # Nn(0,k*a12,Z) = k*a12*diffhankelh1(0,k*a12)*besselj(0,Z) - Z*hankelh1(0,k*a12)*diffbesselj(0,Z)


    Z0 = T(150)

    Qp_eff(Z) = log(Q0(Z/as[1,1],1,1,0,0))/(Z - k_eff*as[1,1])
    Fp_eff = (k_eff*as[1,1] + k*as[1,1]) * exp((T(1.0)/(T(2)*pi*im)) * (
        sum(Fun(Qp_eff,(-Z0)..(as[1,1]*k*(-1.0+1.0im)), num_coefs)) +
        sum(Fun(Qp_eff,(as[1,1]*k*(-1.0+1.0im))..(-as[1,1]*k*(-1.0+1.0im)), num_coefs)) +
        sum(Fun(Qp_eff,(-as[1,1]*k*(-1.0+1.0im))..Z0, num_coefs))
    ))

    Qp_a(Z) = log(Q0(Z/as[1,1],1,1,0,0))/(Z - k*as[1,1])
    Fp_a = (k*as[1,1] + k*as[1,1]) * exp((T(1.0)/(T(2)*pi*im)) * (
        sum(Fun(Qp_a,(-Z0)..(as[1,1]*k*(-1.0+1.0im)), num_coefs)) +
        sum(Fun(Qp_a,(as[1,1]*k*(-1.0+1.0im))..(-as[1,1]*k*(-1.0+1.0im)), num_coefs)) +
        sum(Fun(Qp_a,(-as[1,1]*k*(-1.0+1.0im))..Z0, num_coefs))
    ))
    # has been tested against Mathematica for at least one k_eff
    dsF00_eff = T(2)*k_eff*as[1,1]^2 * Q0(k_eff,1,1,0,0) + ((k_eff*as[1,1])^2 - (k*as[1,1])^2)*dQ0(1,1)

    return T(2)*k*as[1,1]*(cos(θin)/cos(θ_eff))*(Fp_eff/Fp_a)/dsF00_eff
end
