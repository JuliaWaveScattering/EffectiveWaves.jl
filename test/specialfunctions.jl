@testset "Testing special functions" begin
    x = 1.1 + 1.1im
    @test sbesselj(1, x) ≈ sin(x)/x^2 - cos(x)/x
    @test shankelh1(1, x) ≈ - exp(im*x) * (x + im) / (x^2)

    x = 3.1 + 1.1im
    n = 1
    @test 2 * diffsbessel(shankelh1,n,x) ≈ shankelh1(n-1, x) - (shankelh1(n, x) + x*shankelh1(n+1, x))/x

    x = 3.1 - 2.1im
    n = 2
    @test 2 * diffsbessel(shankelh1,n,x) ≈ shankelh1(n-1, x) - (shankelh1(n, x) + x*shankelh1(n+1, x))/x
end
