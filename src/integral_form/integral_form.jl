type Scattering_Amplitudes{T<:Real}
  hankel_order::Int # largest hankel order
  x::Vector{T} # spatial mesh
  A_mat::Matrix{Complex{T}} # a matrix of the scattering amplitudes, size(A_mat) = (length(x), 2hankel_order +1)
end

function Scattering_Amplitudes(x::AbstractVector{T}, A_mat::AbstractMatrix{Complex{T}}) where T<:Number
    Scattering_Amplitudes(Int((size(A_mat,2)-1)/2), x, A_mat)
end

function integrate_B_full(n::Int,X, Y0; Y1 =1000000, θin = 0.0, num_coefs = 10000)
    K(Y) = cos(Y*sin(θin) + n*atan2(Y,X))*hankelh1(n,sqrt(X^2+Y^2))
    # approximate function with Chebyshev polynomial (to high precision) then integrate from Y0 to Y1
    return 2.0*(-1.0)^n*sum(Fun(K,Y0..Y1, num_coefs))
end

# Y0 = sqrt(k^a12^2 - X^2)
function integrate_B(n::Int,X, Y0; θin = 0.0, num_coefs = 10000)
    Y1 = max(2000.0*X, 4000.0) # note Y1 is non-dimensional!
    # assymptotically approximate the integral from Y1 to Inf (tested in integrate_hankels.nb)
    Binf = (1.0+1.0im)*exp(im*Y1*(1.0 - sin(θin)))*
        (1.0 + (-1.0)^n*exp(2.0im*Y1*sin(θin))*(1.0 - sin(θin)) + sin(θin))/(sqrt(pi*Y1)*cos(θin)^2)

    return Binf + integrate_B_full(n, X, Y0; Y1=Y1, θin=θin, num_coefs = num_coefs)
end

# for only whole-correction, this doesn't involve an integral
function integrate_S(n::Int,X; θin = 0.0)
    S = 2.0*(im^Float64(n))*exp(-im*n*θin)*exp(im*X*cos(θin))/cos(θin)
    if X<0 S = conj(S) end
    S
end

function scattering_amplitudes_integral_form(ω::T,medium::Medium{T},specie::Specie{T}; kws...) where T<:Number

    k = ω/medium.c
    (x, (MM_quad,b_mat)) = integral_form(ω, medium, specie;  kws...);

    M = Int( (size(b_mat,2) - 1)/2 )
    J = length(collect(x)) - 1

    len = (J + 1) * (2M + 1)
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    As = MM_mat\b
    As_mat = reshape(As, (J+1, 2M+1))

    return Scattering_Amplitudes(M, collect(x), As_mat)
end

function reflection_coefficient_integrated(ω::T, medium::Medium, specie::Specie,
        amps::Scattering_Amplitudes{T} = scattering_amplitudes_integral_form(ω,medium,specie);
        θin::T = 0.0) where T <: Number

    M = amps.hankel_order
    σ =  trap_scheme(amps.x)
    Z = Array{Complex{Float64}}(2M+1);
    for m = 0:M
        Z[m+M+1] = Zn(ω,specie,medium,m)
        Z[M+1-m] = Z[m+M+1]
    end
    R = T(2)*specie.num_density/(cos(θin)*k^2)*sum( im^T(m)*exp(-im*θin*m)*Z[m+M+1]*amps.A_mat[j,m+M+1]*exp(im*amps.x[j]*cos(θin))*σ[j]
    for m=-M:M, j in eachindex(amps.x))

    return R
end

function integral_form(ω::Float64, medium::Medium, specie::Specie;
        θin::Float64 = 0.0,
        x::AbstractVector = [0.], mesh_points::Int = 201,
        hankel_order = maximum_hankel_order(ω, medium, [specie]; tol=1e-3))

    k = ω/medium.c
    ak = real(k*specie.r);
    M = hankel_order;

    # estimate a large enough mesh
    if x == [0.]
        k_eff = wavenumber_low_volfrac(ω, medium, [specie])
        max_x = 10.0*k/imag(k_eff) # at this A ≈ exp(-10) ≈ 4.5e-5
        J = mesh_points - 1
        h = ak/Int(round(J*ak/max_x));
        x = OffsetArray((0:J)*h, 0:J)
    else
        J = length(collect(x)) - 1
        h = x[2] - x[1]
    end

    Z = OffsetArray{Complex{Float64}}(-M:M);
    for m = 0:M
        Z[m] = Zn(ω,specie,medium,m)
        Z[-m] = Z[m]
    end

    # integration scheme: trapezoidal
    σ =  OffsetArray(trap_scheme(collect(x)), 0:J)
    PQ_quad = intergrand_kernel(x; ak = ak, θin = θin, M = M);

    MM_quad = [
        specie.num_density*Z[n]*σ[j]*PQ_quad[l+1,m+M+1,j+1,n+M+1] + k^2*( (m==n && j==l) ? 1.0+0.0im : 0.0+0.0im)
    for  l=0:J, m=-M:M, j=0:J, n=-M:M];

    b_mat = [ -k^2*exp(im*x[l]*cos(θin))*exp(im*m*(pi/2.0 - θin)) for l = 0:J, m = -M:M]

    return (x, (MM_quad,b_mat))
end

function intergrand_kernel(x::AbstractVector; ak::Float64 = 1.0, θin::Float64 = 0.0,
        M::Int = 2, num_coefs::Int = 10000)

    dx = x[2] - x[1]
    J = length(collect(x)) -1

    if !(typeof(x) <: OffsetArray)
        if J*dx != x[end] warn("Unexpected x = $x.") end
        x = OffsetArray((0:J)*dx, 0:J)
    end
    if !(Int(floor(ak/dx)) ≈ ak/dx)
        warn("There are no mesh points exactly on-top of the intergrands kinks. This could lead to poor accuracy.")
    end
    p = min(Int(floor(ak/dx)),J)
    X = OffsetArray((-J:J)*dx, -J:J)

    B = OffsetArray{Complex{Float64}}(-p:p, -2M:2M);
    for j = -p:p, m = -2M:2M
        if ak^2 -X[j]^2 < -dx^2 error("evaluating B in the wrong domain") end
        B[j,m] = integrate_B(m, X[j], sqrt(abs(ak^2 -X[j]^2)); θin = θin, num_coefs=num_coefs)
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


# plot(collect(x),[real(A_mat[:,M+1]),imag(A_mat[:,M+1])])
# plot(collect(x),abs.(A_mat[:,M+1]))
