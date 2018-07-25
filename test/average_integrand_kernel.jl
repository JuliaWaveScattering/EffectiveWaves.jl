#tests the values of the integrand against results from the Mathematica file integrate_hankels.nb, which used no approximations.

@testset "Compare integrand values against brute force Mathematica result" begin

    h=0.02; max_x=2.0;
    J = Int(round(max_x/h))
    x = (0:J).*h
    i1 = 1 + Int(round(1/h)) # for x1 = 1
    i2 = 1 + Int(round(2/h)) # for x1 = 2

    # both commented tests below should also pass, but am commenting to save time
    # θin=0.; M=0; n=0;

    # # math was produced by Mathematica
    # math00 = [172.31259187840652, 143.34097100449316, 172.31259187840536]
    # intergrand_quad = intergrand_kernel(x; θin=θin, M=M, num_coefs = 2000);
    # julia00 = [sum(abs.(intergrand_quad[i,M+1,:,M+1+n])) for i in [1,i1,i2]];
    # @test maximum(abs.(1.0 .- julia00./math00)) < 1e-6

    # math0p50 = [192.61985415285892, 155.98244729650008, 192.61985415285886]
    # θin=0.5;
    # intergrand_quad = intergrand_kernel(x; θin=θin, M=M, num_coefs = 3000);
    # julia0p50 = [sum(abs.(intergrand_quad[i,M+1,:,M+1+n])) for i in [1,i1,i2]];
    # @test maximum(abs.(1.0 .- julia0p50./math0p50)) < 3e-6

    math0p32 = [241.5144625759003, 211.13367934660897, 181.54133990599098]
    θin=0.3; M=2; n=2;
    intergrand_quad = intergrand_kernel(x; θin=θin, M=M, num_coefs = 2000);
    julia0p32 = [sum(abs.(intergrand_quad[i,M+1,:,M+1+n])) for i in [1,i1,i2]];
    @test maximum(abs.(1.0 .- julia0p32./math0p32)) < 4e-6

end
