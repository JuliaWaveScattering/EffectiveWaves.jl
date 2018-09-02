#  Without using ApproxFun
# function integrate_B_full(n::Int, X, Y0; Y1 =1000000, θin = 0.0, num_coefs = 1000000)
#     Ys = linspace(Y0,Y1,num_coefs);
#     σs = integration_scheme(Ys; scheme = :trapezoidal)
#     K(Y) = cos(Y*sin(θin) + n*atan2(Y,X))*hankelh1(n,sqrt(X^2+Y^2))
#     return 2.0*(-1.0)^n*sum(K.(Ys).*σs)
# end

function integrate_B_full(n::Int, X, Y0; Y1 =1000000, θin = 0.0, num_coefs = 10000)
    K(Y) = cos(Y*sin(θin) + n*atan2(Y,X))*hankelh1(n,sqrt(X^2+Y^2))
    # approximate function with Chebyshev polynomial (to high precision) then integrate from Y0 to Y1
    return 2.0*(-1.0)^n*sum(Fun(K,Y0..Y1, num_coefs))
end

# Y0 = sqrt(k^a12^2 - X^2)
function integrate_B(n::Int, X, Y0; θin = 0.0, num_coefs = 10000)
    Y1 = max(2000.0*X, 4000.0) # note Y1 is non-dimensional!
    # assymptotically approximate the integral from Y1 to Inf (tested in integrate_hankels.nb)
    Binf = (1.0+1.0im)*exp(im*Y1*(1.0 - sin(θin)))*
        (1.0 + (-1.0)^n*exp(2.0im*Y1*sin(θin))*(1.0 - sin(θin)) + sin(θin))/(sqrt(pi*Y1)*cos(θin)^2)

    return Binf + integrate_B_full(n, X, Y0; Y1=Y1, θin=θin, num_coefs = num_coefs)
end

# for only whole-correction, this doesn't involve an integral
function integrate_S(n::Int, X::T; θin::T = 0.0) where T <: AbstractFloat
    S = 2.0*(im^T(n))*exp(-im*n*θin)*exp(im*X*cos(θin))/cos(θin)
    if X<0 S = conj(S) end
    S
end

function BS_matrices(X::AbstractVector{T}, a12k::T; θin::T = 0.0,
        M::Int = 2, num_coefs::Int = 10000) where T<:AbstractFloat

    dX = X[2] - X[1]
    J = length(collect(X)) -1
    q = min(Int(floor(a12k/dX)),J)
    X = OffsetArray((-J:J)*dX, -J:J)

    B = OffsetArray{Complex{Float64}}(-q:q, -2M:2M);
    for j = -q:q, m = -2M:2M
        if a12k^2 - X[j]^2 < - dX^2 error("evaluating B in the wrong domain") end
        B[j,m] = integrate_B(m, X[j], sqrt(abs(a12k^2 -X[j]^2)); θin = θin, num_coefs=num_coefs)
    end
    S = OffsetArray{Complex{Float64}}(-J:J, -2M:2M);
    for j = -J:J, m = -2M:2M
        S[j,m] = integrate_S(m, X[j]; θin = θin)
    end

    return (B,S)
end

function intergrand_kernel(X::AbstractVector{T}, a12k::T; M::Int = 2,
        scheme::Symbol = :trapezoidal, kws...) where T<:AbstractFloat

    dX = X[2] - X[1]
    J = length(collect(X)) -1
    if abs(J*dX - X[end])/X[end] > 1e-10
        error("Expected X to be uniformly spaced from 0. Instead got X= $X.")
    end
    if ( abs(Int(round(a12k/dX)) - a12k/dX) > 1e-10 )
        warn("There are no mesh points exactly on-top of the intergrands kinks. This could lead to poor accuracy.")
    end

    q = min(Int(floor(a12k/dX)),J)
    if q == 0 warn("Mesh element larger than ka12. This is only accurate for low frequency.") end

    B, S = BS_matrices(X, a12k; M = M, kws...)

    # an numerical integration scheme for domain over X
    σ = OffsetArray(integration_scheme(X; scheme=scheme), 0:J)
    # an array of numerical integration schemes for the varying domain |j-l| <= q
    σs = OffsetArray([
        OffsetArray(
            integration_scheme(X[max(l-q,0)+1:min(l+q,J)+1]; scheme=scheme),
        max(l-q,0):min(l+q,J))
    for l=0:J], 0:J)

    intergrand_quad = [
        begin
            P = σ[j] * S[j-l,n-m]
            Q = (abs(j-l) <= q) ? σs[l][j] * (B[j-l,n-m] - S[j-l,n-m]) : zero(Complex{T})
            P + Q
        end
    for l=0:J, m=-M:M, j=0:J, n=-M:M]

    return intergrand_quad
end
