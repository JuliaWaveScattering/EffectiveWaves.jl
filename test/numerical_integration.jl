using EffectiveWaves, Test

@testset "test different numerical integration methods" begin

    using ApproxFun
    import SpecialFunctions: hankelh1

    # function imitates kernal in Fredholm equation that is solved for AverageWave
    n = 1; X = 0.9; θin = 0.2;
    num_coefs = 10000;
    tol = (1.0/num_coefs)^2

    K(Y) = cos(Y*sin(θin) + n*atan(Y,X))*hankelh1(n,sqrt(X^2+Y^2))

    # approximate function with Chebyshev polynomial (to high precision) then integrate from Y0 to Y1
    Y0 = 0.0; Y1 = 1000.0;
    s = sum(Fun(K,Y0..Y1, num_coefs))

    # compare with basic discrete integration
    Ys = LinRange(Y0,Y1,2*num_coefs+1);
    σ = integration_scheme(Ys);
    @test abs(sum(K.(Ys).*σ) - s)/abs(s) < 500*tol # trapezoidal scheme not so accurate

    σ = integration_scheme(Ys; scheme=:simpson);
    @test abs(sum(K.(Ys).*σ) - s)/abs(s) < tol

    integration_scheme(Ys; scheme=:fish);
    @test true
end
