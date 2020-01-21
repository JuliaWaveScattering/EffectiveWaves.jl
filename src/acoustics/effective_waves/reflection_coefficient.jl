"The average reflection coefficient"
function reflection_coefficient(ω::T, wave_eff::EffectiveWave{T}, medium::Acoustic{T,2}, species::Species{T,2};
        θin::T = zero(T), x::T = zero(T), kws...) where T<:Number

    k = ω/medium.c
    θ_ref = pi - wave_eff.θ_eff - θin
    S = length(species)
    ho = wave_eff.basis_order

    kθ = (k*cos(θin) + wave_eff.k_eff*cos(wave_eff.θ_eff))
    R = 2.0im / (k*cos(θin) * kθ)
    R = R*sum(
        exp(im*n*θ_ref + im*x*kθ) * number_density(species[l]) *
        wave_eff.amplitudes[n+ho+1,l]
    for n=-ho:ho, l=1:S)

    return R
end

"The average reflection coefficient"
function wienerhopf_reflection_coefficient(ω::T, medium::Acoustic{T,2}, species::Species{T};
        tol::T = T(1e-7),
        θin::T = zero(T),
        basis_order::Int = 0,
        num_coefs::Int = 20000
    ) where T<:Number

    k = ω/medium.c
    ho = basis_order

    t_vecs = get_t_matrices(medium, species, ω, ho)

    as = [
        outer_radius(s1) * s1.exclusion_distance + outer_radius(s2) * s2.exclusion_distance
    for s1 in species, s2 in species]

    sToS(s,j::Int,l::Int) = (real(s) >= 0) ? sqrt(s^2 + (k*as[j,l]*sin(θin))^2) : -sqrt(s^2 + (k*as[j,l]*sin(θin))^2)

    function Ψ(s,j,l,m,n)
        (s^T(2) - (k*as[j,l]*cos(θin))^T(2)) * (n == m ? T(1) : T(0)) * (j == l ? T(1) : T(0)) +
        T(2) * as[j,l]^T(2) * pi * number_density(species[l]) * t_vecs[l][m+ho+1,m+ho+1] *
        kernelN(n-m,k*as[j,l], sToS(s,j,l))
    end

    # kernelN(0,k*a12,Z) = k*a12*diffhankelh1(0,k*a12)*besselj(0,Z) - Z*hankelh1(0,k*a12)*diffbesselj(0,Z)
    q(s,j,l,m,n) = Ψ(s,j,l,m,n) / (s^T(2) - (k*as[j,l]*cos(θin))^T(2))

    Zs = LinRange(T(100),1/(10*tol),3000)
    maxZ = Zs[findfirst(Z -> abs(log(q(Z,1,1,0,0))) < 10*tol, Zs)]

    function Ψp(s, maxZ::T = maxZ, num_coefs::Int = num_coefs)
        Q(z) = log(q(z,1,1,0,0))/(z - s)
        xp = as[1,1]*k*cos(θin)*(-1.0+1.0im)
        q_pos = exp(
            (T(1.0)/(T(2)*pi*im)) * (
                sum(Fun(Q, Segment(-maxZ,xp), num_coefs)) +
                sum(Fun(Q, Segment(xp,-xp), num_coefs)) +
                sum(Fun(Q, Segment(-xp,maxZ), num_coefs))
            )
        )
        return (s + k*as[1,1]*cos(θin)) * q_pos
    end

    function Ψm(s, maxZ::T = maxZ, num_coefs::Int = num_coefs)
        Q(z) = log(q(z,1,1,0,0))/(z - s)
        xm = as[1,1]*k*cos(θin)*(-1.0+0.5im)
        a1 = T(0)
        q_neg = exp(
            -(T(1.0)/(T(2)*pi*im)) * (
                sum(Fun(Q, Segment(-maxZ,xm+a1), num_coefs)) +
                sum(Fun(Q, Segment(xm+a1,-xm+a1), num_coefs)) +
                sum(Fun(Q, Segment(-xm+a1, maxZ), num_coefs))
            )
        )
        return (s - k*as[1,1]*cos(θin)) * q_neg
    end

    x = -as[1,1]*k*cos(θin)*(-1.0+0.75im)

    # abs(Fp(x,maxZ,num_coefs) - Fp(x,maxZ, Int(round(num_coefs*1.1)))) / abs(Fp(x,maxZ,num_coefs))
    # abs(Fm(x,maxZ,num_coefs) - Fm(x,maxZ, Int(round(num_coefs*1.1)))) / abs(Fm(x,maxZ,num_coefs))

    x2 = as[1,1]*k*cos(θin)*(1.0+2.75im)

    err = abs(Ψp(x,maxZ,num_coefs) * Ψm(x,maxZ,num_coefs) - Ψ(x,1,1,0,0))/abs(Ψ(x,1,1,0,0))
    if err > tol
        @warn "Analytic split recovers original function with $err tolerance, instead of the specified tolernace: $tol"
    end

    R = Ψ(k*as[1,1]*cos(θin),1,1,0,0) / (Ψp(k*as[1,1]*cos(θin),maxZ,num_coefs))^2

    return R
end

# function F0(S,j,l,m,n)
#     (S^T(2) - (k*as[j,l])^T(2)) * (n == m ? T(1) : T(0)) * (j == l ? T(1) : T(0)) +
#     T(2) * as[j,l]^T(2) * pi* number_density(species[l]) *t_vecs[l][m+ho+1] * kernelN(n-m,k*as[j,l],S)
# end

# kernelN(0,k*a12,Z) = k*a12*diffhankelh1(0,k*a12)*besselj(0,Z) - Z*hankelh1(0,k*a12)*diffbesselj(0,Z)
# Q0(S,j,l,m,n) = F0(S,j,l,m,n) / (S^T(2) - (k*as[j,l])^T(2))

# function F0p(S, maxZ::T = maxZ, num_coefs::Int = num_coefs)
#     Q(Z) = log(Q0(Z,1,1,0,0))/(Z - S)
#     xp = as[1,1]*k*(-1.0+1.0im)
#     (S + k*as[1,1]) * exp(
#         (T(1.0)/(T(2)*pi*im)) * (
#             sum(Fun(Q, Segment(-maxZ, xp), num_coefs)) +
#             sum(Fun(Q, Segment(xp,-xp), num_coefs)) +
#             sum(Fun(Q, Segment(-xp,maxZ), num_coefs))
#         )
#     )
# end
