using ApproxFun
using OffsetArrays
using EffectiveWaves

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

function integral_form()

    medium = Medium(1.0,1.0+0.0im)
    specie = Specie(ρ=0.1,r=1., c=0.1, volfrac=0.2)

    k=1.;
    a=specie.r;
    ω = real(k*medium.c)
    k_eff = wavenumber_low_volfrac([ω], medium, [specie])
    # physical parameters
    θin = 0.0
    max_x = 10.0*k/imag(k_eff[1]) # at this A ≈ exp(-10) ≈ 4.5e-5

    # discretization parameters
    M = 3;
    h = a*k/15.0;
    J = Int(round(max_x/h)) # choose an even number for Integration schemes
    x = OffsetArray{Float64}(0:J)
    x[0:J] = (0:J)*h


    Z = OffsetArray{Complex{Float64}}(-M:M);
    for m = 0:M
        Z[m] = Zn(ω,specie,medium,m)
        Z[-m] = Z[m]
    end

    # integration scheme: trapezoidal
    σ =  OffsetArray(trap_scheme(collect(x); xn=max_x), 0:J)

    PQ_quad = intergrand_kernel(x; ak = a*k, θin = θin, M = M);

    MM_quad = [
        specie.num_density*Z[n]*σ[j]*PQ_quad[l+1,m+M+1,j+1,n+M+1] + k^2*( (m==n && j==l) ? 1.0+0.0im : 0.0+0.0im)
    for  l=0:J, m=-M:M, j=0:J, n=-M:M]

    len = (J + 1) * (2M + 1)
    MM_mat = reshape(MM_quad, (len, len))

    b_mat = [ -k^2*exp(im*x[l]*cos(θin))*exp(im*m*(pi/2.0 - θin)) for l = 0:J, m = -M:M]
    b = reshape(b_mat, (len))

    As = MM_mat\b
    # b ≈ [ sum(MM_mat[a,b]*A[b] for b=1:len) for a=1:len]

    A_mat = reshape(As, (J+1, 2M+1))
    # b_mat ≈ [ sum(MM_quad[l,m,j,n]*A_mat[j,n] for j=1:(J+1), n=1:(2M+1)) for l=1:(J+1), m=1:(2M+1)]

    return x, A_mat
end


function intergrand_kernel(x::OffsetVector; ak::Float64 = 1.0, θin::Float64 = 0.0, M::Int = 2, num_coefs::Int = 10000)

    # discretization parameters
    h = x[2] - x[1]
    J = Int(round(x[end]/h))
    if !((0:J).*h ≈ x.parent)  warn("unexpected x = $x") end

    p = min(Int(floor(ak/h)),J)
    X = OffsetArray{Float64}(-J:J)
    X[-J:J] = (-J:J)*h

    B = OffsetArray{Complex{Float64}}(-p:p, -2M:2M);
    for j = -p:p, m = -2M:2M
        if ak^2 -X[j]^2 < -h^2 error("evaluating B in the wrong domain") end
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

# tests the integration scheme
function check_integration(k_eff::Complex{Float64} = 1.0+1.0im; k=1.,a=1., h = a*k/55., max_x = 6.0)

    A(n,x) = im^Float64(n)*exp(-im*n*θin)*exp(im*x*k_eff)
    k=1.; a=1.;
    # physical parameters
    θin = 0.3

    # discretization parameters
    M = 3;
    J = Int(round(max_x/h)) # choose an even number for Integration schemes
    x = OffsetArray{Float64}(0:J)
    x[0:J] = (0:J)*h

    # trapezoidal
    σ =  OffsetArray(trap_scheme(collect(x); xn=max_x), 0:J)

    PQ_quad = intergrand_kernel(x; ak = a*k, θin = θin, M = M)

    # Apply integration scheme
    for j=0:J
        PQ_quad[:,:,j+1,:] *= σ[j]
    end
    # PQ_quad = [σ[j]*PQ_quad[l+1,m+M+1,j+1,n+M+1] for l=0:J, m=-M:M, j=0:J, n=-M:M]

    len = (J + 1) * (2M + 1)
    PQ_mat = reshape(PQ_quad, (len, len))

    A_mat = [ A(m,x[l]) for l = 0:J, m = -M:M]
    As = reshape(A_mat, (len))
    As_integrated = reshape(PQ_mat*As, (J+1,2M+1))

    i1 = Int(round(1/h)) # for x1 = 1
    i2 = Int(round(2/h)) # for x1 = 2
    [As_integrated[1,:], As_integrated[1+i1,:], As_integrated[1+i2,:]]
end


# plot(collect(x),[real(A_mat[:,M+1]),imag(A_mat[:,M+1])])
# plot(collect(x),abs.(A_mat[:,M+1]))
