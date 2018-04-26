using ApproxFun
using OffsetArrays
using EffectiveWaves

function integrate_B_full(n::Int,X, Y0; Y1 =1000000, θin = 0.0)
    K(Y) = cos(Y*sin(θin) + n*atan2(Y,X))*hankelh1(n,sqrt(X^2+Y^2))
    # approximate function with Chebyshev polynomial (to high precision) then integrate from Y0 to Y1
    return 2.0*(-1.0)^n*sum(Fun(K,Y0..Y1))
end

# for only whole-correction, this doesn't involve an integral
function integrate_S(n::Int,X; θin = 0.0)
    2.0*(im^Float64(n))*exp(-im*n*θin)*exp(im*X*cos(θin))/cos(θin)
end

# Y0 = sqrt(k^a12^2 - X^2)
function integrate_B(n::Int,X, Y0; θin = 0.0)
    Y1 = max(2000.0*X, 4000.0) # note Y1 is non-dimensional!
    # assymptotically approximate the integral from Y1 to Inf (tested in integrate_hankels.nb)
    Binf = (1.0+1.0im)*exp(im*Y1*(1.0 - sin(θin)))*
        (1.0 + (-1.0)^n*exp(2.0im*Y1*sin(θin))*(1.0 - sin(θin)) + sin(θin))/(sqrt(pi*Y1)*cos(θin)^2)

    return Binf + integrate_B_full(n, X, Y0; Y1=Y1, θin=θin)
end


# "pre-calculate a matrix of Zn's"
# function Zn_matrix(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; hankel_order = 3) where T <: Number
#     Zs = OffsetArray{Complex{Float64}}(1:length(species), -hankel_order:hankel_order)
#     for i = 1:length(species), n = 0:hankel_order
#         Zs[i,n] = Zn(ω,species[i],medium,n)
#         Zs[i,-n] = Zs[i,n]
#     end
#     return Zs
# end

function integral_form(ks)

    medium = Medium(1.0,1.0+0.0im)
    specie = Specie(ρ=0.5,r=0.4, c=0.8, volfrac=0.2)

    M = 3;
    J = 200; # choose an even number for Integration schemes
    k=1.0;
    ω = real(k*medium.c)
    a=0.2;
    θin = 0.0
    h = a*k/10.0;

    p = Int(floor(a*k/h))
    x = OffsetArray{Float64}(0:J)
    x[0:J] = (0:J)*h
    X = OffsetArray{Float64}(-J:J)
    X[-J:J] = (-J:J)*h

    σ = OffsetArray{Complex{Float64}}(0:J)
    # Simpson's integration
    for j = 0:J
        # σ[j] = isodd(j) ? 4.0 : 2.0
        σ[j] = 1.0
    end
    # σ[0] = σ[J] = 1.0
    # σ = σ*h/3.0
    σ = σ*h

    Z = OffsetArray{Complex{Float64}}(-M:M)
    for m = 0:M
        Z[m] = Zn(ω,specie,medium,m)
        Z[-m] = Z[m]
    end

    B = OffsetArray{Complex{Float64}}(-p:p, -2M:2M)
    for j = -p:p, m = -2M:2M
        B[j,m] = integrate_B(m, X[j], sqrt(a^2*k^2 -X[j]^2); θin = θin)
    end
    S = OffsetArray{Complex{Float64}}(-p:J, -2M:2M)
    for j = -p:J, m = -2M:2M
        S[j,m] = integrate_S(m, X[j]; θin = θin)
    end

    function MM(l,j,m,n)
        P = exp(-im*x[l]*cos(θin))*σ[j]*S[j,n-m]
        Q = (abs(j-l)<=p) ? σ[j]*(B[j-l,n-m] - S[j-l,n-m]) : 0.0+0.0im
        specie.num_density*Z[n]*(P + Q) + k^2*( (m==n && j==l) ? 1.0+0.0im : 0.0+0.0im)
    end

    MM_mat = [MM(l,j,m,n) for  j=0:J, m=-M:M, l=0:J, n=-M:M]
    len = (J + 1) * (2M + 1)
    MM_mat = reshape(MM_mat, (len, len))

    b = [ -k^2*exp(im*X[j]*cos(θin)*exp(im*m*(pi/2.0 - θin))) for j = 0:J, m = -M:M]
    b = reshape(b, (len))

    A = MM_mat\b
    A = reshape(A, (J + 1, 2M+1))

    using Plots; pyplot()
    plot(collect(x),[real(A[:,M+1]),imag(A[:,M+1])])
    plot(collect(x),abs.(A[:,M+1]))
end
