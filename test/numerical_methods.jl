# quick test of numerical methods such as integration schemes and function approximations

using Test, LinearAlgebra, Statistics

@testset "Tests for numerical integration" begin
    n = 100
    x = sort(rand(n))
    ws = trapezoidal_scheme(x)

    # test integral of x 
    @test sum(x .* ws) ≈ (x[end]^2 - x[1]^2)/2

    # test integral of x^2. Has quadratic error
    @test abs(sum((x.^2) .* ws) - (x[end]^3 - x[1]^3)/3) < 4 * (1 / n)^2
end
