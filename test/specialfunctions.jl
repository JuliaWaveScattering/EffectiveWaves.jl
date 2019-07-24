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

@testset "Associated legendre functions" begin
    x = rand(1)[1] * 0.99
    l_max = 2

    # a few associated legendre functions without the Condon-Shortly factor
    Plm_arr = [1,x,sqrt(1-x^2),(3x^2-1)/2, 3x*sqrt(1-x^2),3*(1-x^2)]

    ind_max = Int(sf_legendre_array_index(l_max,l_max+1)) # max index, but why do I need to add 1?
    lm_indices = vcat(associated_legendre_indices(l_max)...)

    @test length(lm_indices) == ind_max

    @test sf_legendre_array(GSL_SF_LEGENDRE_NONE, l_max, x)[1:ind_max] ≈ Plm_arr

    sph_factors = map(lm_indices) do lm
       l = lm[1]
       m = lm[2]
       (-1)^m * sqrt((2l + 1)/(4pi) * factorial(l-m) / factorial(l+m))
    end

    condon_phase = [(-1)^lm[2] for lm in lm_indices]

    @test condon_phase .* sf_legendre_array(GSL_SF_LEGENDRE_SPHARM, l_max, x)[1:ind_max] ≈ sph_factors .* Plm_arr

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
