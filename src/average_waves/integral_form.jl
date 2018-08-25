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

function intergrand_kernel(X::AbstractVector{T}, a12k::T; θin::T = 0.0,
        M::Int = 2, num_coefs::Int = 10000) where T<:AbstractFloat

    dX = X[2] - X[1]
    J = length(collect(X)) -1

    if !(typeof(X) <: OffsetArray)
        if abs(J*dX - X[end])/X[end] > 1e-10  warn("Unexpected X = $X.") end
        X = OffsetArray((0:J)*dX, 0:J)
    end
    if ( abs(Int(floor(a12k/dX)) - a12k/dX) > 1e-5 )
        warn("There are no mesh points exactly on-top of the intergrands kinks. This could lead to poor accuracy.")
    end
    p = min(Int(floor(a12k/dX)),J)
    X = OffsetArray((-J:J)*dX, -J:J)

    B = OffsetArray{Complex{Float64}}(-p:p, -2M:2M);
    for j = -p:p, m = -2M:2M
        if a12k^2 - X[j]^2 < - dX^2 error("evaluating B in the wrong domain") end
        B[j,m] = integrate_B(m, X[j], sqrt(abs(a12k^2 -X[j]^2)); θin = θin, num_coefs=num_coefs)
    end
    S = OffsetArray{Complex{Float64}}(-J:J, -2M:2M);
    for j = -J:J, m = -2M:2M
        S[j,m] = integrate_S(m, X[j]; θin = θin)
    end
    function intergrand(l,j,m,n)
        P = S[j-l,n-m]
        Q = (abs(j-l)<= p) ? (B[j-l,n-m] - S[j-l,n-m]) : 0.0+0.0im
        P + Q
    end

    intergrand_quad = [intergrand(l,j,m,n) for  l=0:J, m=-M:M, j=0:J, n=-M:M]

    return intergrand_quad
end
