@testset "spherical bessel functions" begin
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

@testset "Gaunt and Wigner symbols" begin
    l1 = rand(1:100)
    l2 = rand(1:100)
    l3 = rand(abs(l1-l2):2:(l1+l2)) # guarantees that iseven(l1+l2+l3)

    @test iseven(l1+l2+l3)

    m1 = rand(-l1:l1)
    m2 = rand(-l2:l2)
    m3 = m1-m2

    @test_throws(DomainError,gaunt_coefficients(l1,2*l1,l2,m2,l3,m3))
    @test_throws(DomainError,gaunt_coefficients(l1,m1,l2,m2,0.1,m3))
end
