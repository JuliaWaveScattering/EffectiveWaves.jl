"The average reflection coefficient"
function reflection_coefficient(ω::T, wave_eff::EffectivePlaneWaveMode{T}, psource::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{2,Halfspace{T,2}};
        x::T = zero(T), kws...) where T<:Number

    k = ω / psource.medium.c
    θ_eff = transmission_angle(wave_eff,material)
    θin = transmission_angle(psource,material)

    θ_ref = pi - θ_eff - θin
    S = length(material.species)
    ho = wave_eff.basis_order

    kcos_eff = wave_eff.wavenumber * dot(- conj(material.shape.normal), wave_eff.direction)
    kcos_in = k * dot(- conj(material.shape.normal), psource.direction)

    kθ = kcos_in + kcos_eff
    R = 2.0im / (kcos_in * kθ)
    R = R*sum(
        exp(im*n*θ_ref + im*x*kθ) * number_density(material.species[l]) *
        wave_eff.amplitudes[n+ho+1,l]
    for n=-ho:ho, l=1:S)

    return R
end

"The average reflection coefficient"
function wienerhopf_reflection_coefficient(ω::T, psource::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{2,Halfspace{T,2}};
        tol::T = T(1e-7),
        basis_order::Int = 0,
        num_coefs::Int = 20000
    ) where T<:Number

    k = ω/psource.medium.c
    ho = basis_order

    t_vecs = get_t_matrices(psource.medium, material.species, ω, ho)

    θin = transmission_angle(psource, material)
    kcos = k*cos(θin)
    ksin = k*sin(θin)

    as = [
        outer_radius(s1) * s1.exclusion_distance + outer_radius(s2) * s2.exclusion_distance
    for s1 in material.species, s2 in material.species]

    sToS(s,j::Int,l::Int) = (real(s) >= 0) ? sqrt(s^2 + (as[j,l]*ksin)^2) : -sqrt(s^2 + (as[j,l]*ksin)^2)

    function Ψ(s,j,l,m,n)
        (s^T(2) - (as[j,l]*kcos)^T(2)) * (n == m ? T(1) : T(0)) * (j == l ? T(1) : T(0)) +
        T(2) * as[j,l]^T(2) * pi * number_density(material.species[l]) * t_vecs[l][m+ho+1,m+ho+1] *
        kernelN2D(n-m,k*as[j,l], sToS(s,j,l))
    end

    # kernelN2D(0,k*a12,Z) = k*a12*diffhankelh1(0,k*a12)*besselj(0,Z) - Z*hankelh1(0,k*a12)*diffbesselj(0,Z)
    q(s,j,l,m,n) = Ψ(s,j,l,m,n) / (s^T(2) - (as[j,l]*kcos)^T(2))

    Zs = LinRange(T(100),1/(10*tol),3000)
    maxZ = Zs[findfirst(Z -> abs(log(q(Z,1,1,0,0))) < 10*tol, Zs)]

    function Ψp(s, maxZ::T = maxZ, num_coefs::Int = num_coefs)
        Q(z) = log(q(z,1,1,0,0))/(z - s)
        xp = as[1,1]*kcos * (-1.0+1.0im)
        q_pos = exp(
            (T(1.0)/(T(2)*pi*im)) * (
                sum(Fun(Q, Segment(-maxZ,xp), num_coefs)) +
                sum(Fun(Q, Segment(xp,-xp), num_coefs)) +
                sum(Fun(Q, Segment(-xp,maxZ), num_coefs))
            )
        )
        return (s + as[1,1] * kcos) * q_pos
    end

    function Ψm(s, maxZ::T = maxZ, num_coefs::Int = num_coefs)
        Q(z) = log(q(z,1,1,0,0))/(z - s)
        xm = as[1,1] * kcos * (-1.0+0.5im)
        a1 = T(0)
        q_neg = exp(
            -(T(1.0)/(T(2)*pi*im)) * (
                sum(Fun(Q, Segment(-maxZ,xm+a1), num_coefs)) +
                sum(Fun(Q, Segment(xm+a1,-xm+a1), num_coefs)) +
                sum(Fun(Q, Segment(-xm+a1, maxZ), num_coefs))
            )
        )
        return (s - as[1,1] * kcos) * q_neg
    end

    x = -as[1,1] * kcos * (-1.0+0.75im)

    # abs(Fp(x,maxZ,num_coefs) - Fp(x,maxZ, Int(round(num_coefs*1.1)))) / abs(Fp(x,maxZ,num_coefs))
    # abs(Fm(x,maxZ,num_coefs) - Fm(x,maxZ, Int(round(num_coefs*1.1)))) / abs(Fm(x,maxZ,num_coefs))

    x2 = as[1,1] * kcos * (1.0+2.75im)

    err = abs(Ψp(x,maxZ,num_coefs) * Ψm(x,maxZ,num_coefs) - Ψ(x,1,1,0,0))/abs(Ψ(x,1,1,0,0))
    if err > tol
        @warn "Analytic split recovers original function with $err tolerance, instead of the specified tolernace: $tol"
    end

    R = field(psource,zeros(T,2),ω) * Ψ(as[1,1] * kcos,1,1,0,0) / (Ψp(as[1,1] * kcos,maxZ,num_coefs))^2

    return R
end

# function F0(S,j,l,m,n)
#     (S^T(2) - (k*as[j,l])^T(2)) * (n == m ? T(1) : T(0)) * (j == l ? T(1) : T(0)) +
#     T(2) * as[j,l]^T(2) * pi* number_density(species[l]) *t_vecs[l][m+ho+1] * kernelN2D(n-m,k*as[j,l],S)
# end

# kernelN2D(0,k*a12,Z) = k*a12*diffhankelh1(0,k*a12)*besselj(0,Z) - Z*hankelh1(0,k*a12)*diffbesselj(0,Z)
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
