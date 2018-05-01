"tests the values of the integrand against results from the Mathematica file integrate_hankels.nb, which used no approximations."
function test_integrand()

    h=0.02; max_x=2.0;
    J = Int(round(max_x/h))
    x = (0:J).*h
    i1 = 1 + Int(round(1/h)) # for x1 = 1
    i2 = 1 + Int(round(2/h)) # for x1 = 2

    θin=0.; M=0; n=0;

    # math was produced by Mathematica
    math00 = [172.31259187840652, 143.34097100449316, 172.31259187840536]
    intergrand_quad = intergrand_kernel(x; θin=θin, M=M, num_coefs = 5000);
    julia00 = [sum(abs.(intergrand_quad[i,M+1,:,M+1+n])) for i in [1,i1,i2]];
    @test maximum(abs.(1.0 .- julia00./math00)) < 1e-6


    math0p50 = [192.61985415285892, 155.98244729650008, 192.61985415285886]
    θin=0.5;
    intergrand_quad = intergrand_kernel(x; θin=θin, M=M, num_coefs = 10000);
    julia0p50 = [sum(abs.(intergrand_quad[i,M+1,:,M+1+n])) for i in [1,i1,i2]];
    @test maximum(abs.(1.0 .- julia0p50./math0p50)) < 3e-6

    math0p32 = [241.5144625759003, 211.13367934660897, 181.54133990599098]
    θin=0.3; M=2; n=2;
    intergrand_quad = intergrand_kernel(x; θin=θin, M=M, num_coefs = 5000);
    julia0p32 = [sum(abs.(intergrand_quad[i,M+1,:,M+1+n])) for i in [1,i1,i2]];
    @test maximum(abs.(1.0 .- julia0p32./math0p32)) < 4e-6

end

# ints = [ check_integration(1.0+1.0im; h = 1./n) for n=10:30:310]

function test_check_integration()
#from Mathematics, had several errors when running
math = [
        [65.0544 - 146.891im, 24.1498 - 8.37225im, -5.10012 + 3.74221im
            , 0.940428 - 3.4201im, 7.25312 - 0.0534466im, 22.9058 + 38.434im, -49.5574 + 92.4377im]
        ,[-32.7224 + 1.52201im, -12.033 - 16.2499im, 5.93884 - 10.8663im
            , 8.51838 + 3.68509im, -4.94285 + 9.32462im, -12.0994 - 5.48456im, 3.6755 - 13.9906im]
        ,[-9.82076 - 18.6418im, 10.2188 - 12.2782im, 11.8699 + 5.34705im
            , -2.698 + 12.0431im, -13.1158 - 0.218682im, -2.92199 - 13.8296im, 12.9993 - 6.4169im]
     ]
#from julia
using JLD
using Plots; pyplot()
ints =  first(values(load("integrated_As.jld")))
data = [ [ints[i][j][k] for i in eachindex(ints)] for j=1:3, k=1:7];
# i=2;j=3;k=5; ints[i][j][k] == data[j,k][i]
# x = 0. => j = 1
M = 3;

hs =  [ 1./n for n=10:30:310]
j = 3
 plot()
 for k=1:7
     plot!(hs[4:end], abs.(data[j,k][4:end]), label = "m = $(k - M - 1)")
     scatter!([hs[end]], [abs(math[j][k])])
 end
  gui()

 plot()
 for k=1:7
     plot!(hs[4:end], 1 - abs.(data[j,k][4:end])./abs(data[j,k][end]), label = "m = $(k - M - 1)")
 end
  gui()
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
