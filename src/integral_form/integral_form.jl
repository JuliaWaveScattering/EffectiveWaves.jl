using ApproxFun
using OffsetArrays
using EffectiveWaves
using Plots; pyplot()

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
    specie = Specie(ρ=0.1,r=0.4, c=0.1, volfrac=0.2)

    k=2.5;
    a=specie.r;
    ω = real(k*medium.c)
    k_eff = wavenumber_low_volfrac([ω], medium, [specie])
    # physical parameters
    θin = 0.0
    max_x = 10.0*k/imag(k_eff[1]) # at this A ≈ exp(-10) ≈ 4.5e-5

    # num_coefs = 2000
    # B_funs = OffsetArray{ApproxFun.Fun}(-M:M)
    # for m = -M:M
    #     B_funs[m] = Fun(
    #         X -> integrate_B(m, X, sqrt(abs(a^2*k^2 -X^2)); θin = θin, num_coefs = num_coefs)
    #     ,0.0..max_x, num_coefs)
    # end

    # discretization parameters
    M = 3;
    h = a*k/15.0;
    J = Int(round(max_x/h)) # choose an even number for Integration schemes

    p = min(Int(floor(a*k/h)),J)
    x = OffsetArray{Float64}(0:J)
    x[0:J] = (0:J)*h
    X = OffsetArray{Float64}(-J:J)
    X[-J:J] = (-J:J)*h

    σ = OffsetArray{Complex{Float64}}(0:J)

    for j = 0:J
        # σ[j] = isodd(j) ? 4.0 : 2.0 # Simpson's integration
        σ[j] = 1.0 # trapezoidal is better suited to discontinuous curves
    end
    σ[0] = σ[J] = 0.5
    # σ = σ*h/3.0
    σ = σ*h

    Z = OffsetArray{Complex{Float64}}(-M:M);
    for m = 0:M
        Z[m] = Zn(ω,specie,medium,m)
        Z[-m] = Z[m]
    end

    B = OffsetArray{Complex{Float64}}(-p:p, -2M:2M);
    for j = -p:p, m = -2M:2M
        if a^2*k^2 -X[j]^2 < -h^2 error("evaluating B in the wrong domain") end
        B[j,m] = integrate_B(m, X[j], sqrt(abs(a^2*k^2 -X[j]^2)); θin = θin, num_coefs = 10000)
    end
    S = OffsetArray{Complex{Float64}}(-p:J, -2M:2M);
    for j = -p:J, m = -2M:2M
        S[j,m] = integrate_S(m, X[j]; θin = θin)
    end

    function MM(l,j,m,n)
        P = exp(-im*x[l]*cos(θin))*σ[j]*S[j,n-m]
        Q = (abs(j-l)<=p) ? σ[j]*(B[j-l,n-m] - S[j-l,n-m]) : 0.0+0.0im
        specie.num_density*Z[n]*(P + Q) + k^2*( (m==n && j==l) ? 1.0+0.0im : 0.0+0.0im)
    end

    MM_quad = [MM(l,j,m,n) for  l=0:J, m=-M:M, j=0:J, n=-M:M]
    #NOTE: sum([MM(l,j,m,n) - MM_quad[l+1,m+M+1,j+1,n+M+1] for  j=0:J, n=-M:M, l=0:J, m=-M:M])


    len = (J + 1) * (2M + 1)
    MM_mat = reshape(MM_quad, (len, len))

    b_mat = [ -k^2*exp(im*X[l]*cos(θin))*exp(im*m*(pi/2.0 - θin)) for l = 0:J, m = -M:M]
    b = reshape(b_mat, (len))

    A = MM_mat\b
    # b ≈ [ sum(MM_mat[a,b]*A[b] for b=1:len) for a=1:len]

    A_mat = reshape(A, (J+1, 2M+1))
    # b_mat ≈ [ sum(MM_quad[l,m,j,n]*A_mat[j,n] for j=1:(J+1), n=1:(2M+1)) for l=1:(J+1), m=1:(2M+1)]

    return x, A_mat
end

# ints = [ check_integration(1.0+1.0im; h = 1./n) for n=10:30:310]


function test_check_integration()
#from Mathematics, had several errors when running
# math = [[63.3587 - 146.953im, 23.5093 - 9.523im, -6.79233 + 3.62019im,
#      0.337901 - 4.59524im, 5.48931 - 0.100657im, 22.2992 + 37.1681im,
#      -51.2365 + 92.3946im], [-34.4328 + 1.42769im, -12.6768 - 17.4311im,
#      4.25781 - 11.0139im, 7.912 + 2.52933im, -6.65733 + 9.23779im,
#      -12.6634 - 6.6802im, 1.94403 - 14.0034im],
#     [-11.4586 - 18.8124im, 9.61476 - 13.4174im, 10.16 + 5.25529im,
#      -3.34734 + 10.8727im, -14.7803 - 0.335775im, -3.52956 - 14.9673im,
#      11.2167 - 6.42339im]]
#from julia
ints = load("integrated_As.jld")

end


# tests the integration scheme
function check_integration(k_eff::Complex{Float64} = 1.0+1.0im; k=1.,a=1., h = a*k/55.)

    k=1.; a=1.;
    # physical parameters
    θin = 0.3
    max_x = 2.0

    # discretization parameters
    M = 3;

    J = Int(round(max_x/h)) # choose an even number for Integration schemes

    p = min(Int(floor(a*k/h)),J)
    x = OffsetArray{Float64}(0:J)
    x[0:J] = (0:J)*h
    X = OffsetArray{Float64}(-J:J)
    X[-J:J] = (-J:J)*h

    σ = OffsetArray{Complex{Float64}}(0:J)
    # trapezoidal
    for j = 0:J
        σ[j] = 1.0
    end
    σ[0] = σ[J] = 0.5
    σ = σ*h

    B = OffsetArray{Complex{Float64}}(-p:p, -2M:2M);
    for j = -p:p, m = -2M:2M
        if a^2*k^2 -X[j]^2 < -h^2 error("evaluating B in the wrong domain") end
        B[j,m] = integrate_B(m, X[j], sqrt(abs(a^2*k^2 -X[j]^2)); θin = θin, num_coefs = 20000)
    end
    S = OffsetArray{Complex{Float64}}(-J:J, -2M:2M);
    for j = -J:J, m = -2M:2M
        S[j,m] = integrate_S(m, X[j]; θin = θin)
    end

    function PQ_integral(l,j,m,n)
        P = S[j-l,n-m]
        Q = (abs(j-l)<= p) ? (B[j-l,n-m] - S[j-l,n-m]) : 0.0+0.0im
        σ[j]*(P + Q)
    end
    PQ_quad = [PQ_integral(l,j,m,n) for  l=0:J, m=-M:M, j=0:J, n=-M:M]
    #NOTE: sum([MM(l,j,m,n) - MM_quad[l+1,m+M+1,j+1,n+M+1] for  j=0:J, n=-M:M, l=0:J, m=-M:M])

    len = (J + 1) * (2M + 1)
    PQ_mat = reshape(PQ_quad, (len, len))

    A(n,x) = im^Float64(n)*exp(-im*n*θin)*exp(im*x*k_eff)
    A_mat = [ A(m,X[l]) for l = 0:J, m = -M:M]
    As = reshape(A_mat, (len))

    As_integrated = reshape(PQ_mat*As, (J+1,2M+1))

    i1 = Int(round(1/h)) # for x1 = 1
    i2 = Int(round(2/h)) # for x1 = 2
    [As_integrated[1,:], As_integrated[1+i1,:], As_integrated[1+i2,:]]
end


# plot(collect(x),[real(A_mat[:,M+1]),imag(A_mat[:,M+1])])
# plot(collect(x),abs.(A_mat[:,M+1]))
