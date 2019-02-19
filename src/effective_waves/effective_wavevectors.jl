
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
#
# "calculate effective transmission angle"
# function transmission_angle(k::Complex{T}, k_eff::Complex{T}, θin; tol = 1e-8) where T<:Number
#     snell(θ::Array{T}) = norm(k*sin(θin) - k_eff*sin(θ[1] + im*θ[2]))
#     result = optimize(snell, [θin,0.]; g_tol= tol^2.0)
#     θ_eff = result.minimizer[1] + im*result.minimizer[2]
# end

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

function wienerhopf_wavevectors(ω::T, k_eff::Complex{T}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:AbstractFloat
    return wienerhopf_wavevectors(ω, [k_eff], medium, species; kws...)
end


"The average effective transmitted wavevectors according to the Wiener-Hopf method.
The function returns an array A, where
AA(x,y,0,1) = A[1,1]*exp(im*k_eff*(cos(θ_eff)*x + sin(θin)*y))
where (x,y) are coordinates in the halfspace  and AA is the ensemble average scattering coefficient. Method currently only implemented for 1 species and for monopole scatterers."
function wienerhopf_wavevectors(ω::T, k_effs::Vector{Complex{T}}, medium::Medium{T}, species::Vector{Specie{T}};
        tol::T = 1e-6, θin::T = 0.0,
        radius_multiplier::T = 1.005,
        hankel_order::Int = 0,
        num_coefs::Int = 10000,
        maxZ::T = T(100)*radius_multiplier*maximum(s.r for s in species) + T(100),
        kws...
    ) where T<:AbstractFloat

    k = ω/medium.c
    ho = hankel_order

    t_vecs = t_vectors(ω, medium, species; hankel_order = ho)
    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]

    # differentiate in S ( = a[j,l] k_eff ) and evaluate at S
    function dSQ0_eff(S,j,l)
        m = 0

        T(2) * as[j,l]^T(2) * pi*species[l].num_density*t_vecs[l][m+ho+1]*(-T(2)*S)/(S^T(2) - (k*as[j,l])^T(2))^2 *
        Nn(0,k*as[j,l],S) +
        T(2) * pi * as[j,l]^T(2) * species[l].num_density*t_vecs[l][m+ho+1]/(S^T(2) - (k*as[j,l])^T(2)) * (
            k*as[j,l]*diffhankelh1(0,k*as[j,l])*diffbesselj(0,S) -
            hankelh1(0,k*as[j,l])*diffbesselj(0,S) -
            S*hankelh1(0,k*as[j,l])*diffbesselj(0,S,2)
        )
    end

    # Nn(0,k*a12,Z) = k*a12*diffhankelh1(0,k*a12)*besselj(0,Z) - Z*hankelh1(0,k*a12)*diffbesselj(0,Z)
    # DZNn(0,k*a12,Z) = k*a12*diffhankelh1(0,k*a12)*diffbesselj(0,Z) - Z*hankelh1(0,k*a12)*diffbesselj(0,Z,2) -

    sToS(s,j::Int,l::Int) = (real(s) >= 0) ? sqrt(s^2 + (k*as[j,l]*sin(θin))^2) : -sqrt(s^2 + (k*as[j,l]*sin(θin))^2)

    function q(s,j,l,m,n)
        (n == m ? T(1) : T(0)) * (j == l ? T(1) : T(0)) +
        T(2) * pi * as[j,l]^T(2) * species[l].num_density * t_vecs[l][m+ho+1] *
        Nn(n-m, k*as[j,l], sToS(s,j,l)) / (s^T(2) - (k*as[j,l]*cos(θin))^T(2))
    end

    Zs = LinRange(T(100),1/(10*tol),3000)
    maxZ = Zs[findfirst(Z -> abs(log(q(Z,1,1,0,0))) < 10*tol, Zs)]

    function Fp(s, maxZ::T = maxZ, num_coefs::Int = num_coefs)
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


    Fp_a = Fp(k*as[1,1]*cos(θin))
    # has been tested against Mathematica for at least one k_eff
    dSF00(S) = (S^2 - (k*as[1,1])^2)*dSQ0_eff(S,1,1) # + T(2) * S * Q0(S,1,1,0,0)
    # last term left out becuase Q0(k_eff,1,1,0,0) = 0.

    return map(k_effs) do k_eff
        θ_eff = transmission_angle(k, k_eff, θin; tol = tol)
        Fp_eff = Fp(k_eff*as[1,1]*cos(θ_eff))
        T(2)*k*as[1,1]*(cos(θin)/cos(θ_eff))*(Fp_eff/Fp_a)/dSF00(k_eff*as[1,1])
    end
end
