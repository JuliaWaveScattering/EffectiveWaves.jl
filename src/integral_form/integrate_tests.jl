# tests the values of the integrand against results from the Mathematica file integrate_hankels.nb, which used no approximations.
function test_integrad()

    # math was produced by Mathematica
    math00 = [172.31259187840652, 143.34097100449316, 172.31259187840536]

    h=0.02; max_x=2.0;
    i1 = Int(round(1/h)) # for x1 = 1
    i2 = Int(round(2/h)) # for x1 = 2

    θin=0.; M=0; n=0;
    (x, integrand_quad) = check_integrand(θin,M; max_x = max_x, h = h)
    julia00 = [sum(abs.(integrand_quad[i,M+1,:,M+1+n])) for i in [1,1+i1,1+i2]]
    maximum(1.0 .- julia00./math00) < 1e-6


    math0p50 = [192.61985415285892, 155.98244729650008, 192.61985415285886]
    θin=0.5;
    (x, integrand_quad) = check_integrand(θin,M; max_x = max_x, h = h)
    julia0p50 = [sum(abs.(integrand_quad[i,M+1,:,M+1+n])) for i in [1,1+i1,1+i2]]
    maximum(1.0 .- julia0p50./math0p50) < 3e-6

    math0p32 = [241.5144625759003, 211.13367934660897, 181.54133990599098]
    θin=0.3; M=2; n=2;
    (x, integrand_quad) = check_integrand(θin,M; max_x = max_x, h = h)
    julia0p32 = [sum(abs.(integrand_quad[i,M+1,:,M+1+n])) for i in [1,1+i1,1+i2]]
    maximum(1.0 .- julia0p32./math0p32) < 3e-6

end


# calculates the values of the integrand
function check_integrand(θin::Float64,M::Int; max_x = 2.0, h = 0.02)
    k=1.0; a=1.0;

    # discretization parameters
    J = Int(round(max_x/h)) # choose an even number for Integration schemes

    p = min(Int(floor(a*k/h)),J)
    x = OffsetArray{Float64}(0:J)
    x[0:J] = (0:J)*h
    X = OffsetArray{Float64}(-J:J)
    X[-J:J] = (-J:J)*h

    B = OffsetArray{Complex{Float64}}(-p:p, -2M:2M);
    for j = -p:p, m = -2M:2M
        if a^2*k^2 -X[j]^2 < -h^2 error("evaluating B in the wrong domain") end
        B[j,m] = integrate_B(m, X[j], sqrt(abs(a^2*k^2 -X[j]^2)); θin = θin, num_coefs = 20000)
    end
    S = OffsetArray{Complex{Float64}}(-J:J, -2M:2M);
    for j = -J:J, m = -2M:2M
        S[j,m] = integrate_S(m, X[j]; θin = θin)
    end
    function PQ_integrand(l,j,m,n)
        P = S[j-l,n-m]
        Q = (abs(j-l)<= p) ? (B[j-l,n-m] - S[j-l,n-m]) : 0.0+0.0im
        P + Q
    end

    PQ_quad = [PQ_integrand(l,j,m,n) for  l=0:J, m=-M:M, j=0:J, n=-M:M]

    return (x, PQ_quad)
end
