using GSL
using LinearAlgebra

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

    lm_indices = associated_legendre_positive_indices(l_max)

    @test sf_legendre_array(GSL_SF_LEGENDRE_NONE, l_max, x)[1:ind_max] ≈ Plm_arr

    sph_factors = map(lm_indices) do lm
       l = lm[1]
       m = lm[2]
       (-1)^m * sqrt((2l + 1)/(4pi) * factorial(l-m) / factorial(l+m))
    end

    condon_phase = [(-1)^lm[2] for lm in lm_indices]

    @test condon_phase .* sf_legendre_array(GSL_SF_LEGENDRE_SPHARM, l_max, x)[1:ind_max] ≈ sph_factors .* Plm_arr

end

@testset "Spherical harmonics" begin
    θ = rand(1)[1] * 0.99
    φ = rand(1)[1] * 0.99

    l_max = 8 # small l_max due to factorial formula below

    ls, ms = spherical_harmonics_indices(l_max)
    sphs = spherical_harmonics(l_max, θ, φ)

    # check special case l == abs(m)
    inds = findall(ls .== abs.(ms))

    for i in inds
        @test sphs[i] ≈ (sign(ms[i]))^ls[i] / (2^ls[i] * factorial(ls[i])) *
            sqrt(factorial(2*ls[i] + 1) / (4pi)) * sin(θ)^ls[i] * exp(im * ms[i] * φ)
    end

    l_max = 30
    ls, ms = spherical_harmonics_indices(l_max)
    sphs = spherical_harmonics(l_max, θ, φ)

    # special case m == 0, reduce to just Legendre polynomials
    Ps = sf_legendre_Pl_array(l_max, cos(θ))
    inds = findall(ms .== 0)

    for l in 0:l_max
        @test sphs[inds][l+1] ≈ sqrt((2l+1)/(4pi)) * Ps[l+1]
    end

    #sphs[inds] .≈ sqrt.((2 .* (0:l_max) .+ 1) ./ (4pi)) .* Ps

    # special case, north pole
    θ = 0.0
    sphs = spherical_harmonics(l_max, θ, φ)

    for i in eachindex(sphs)
        if ms[i] == 0
            @test sphs[i] ≈ sqrt((2*ls[i] + 1)/(4pi))
        else
            @test sphs[i] ≈ 0.0
        end
        i += 1
    end

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
    @test_throws(MethodError,gaunt_coefficients(l1,m1,l2,m2,0.1,m3))

    # the spherical harmonics linearisation formula
    θ, φ = rand(2) .* 0.99

    l_small = 6
    l_max = 2*l_small # needs to be larger than l_small

    ls, ms = spherical_harmonics_indices(l_max);
    Ys = spherical_harmonics(l_max, θ, φ);

    for n2 in 1:(l_small + 1)^2, n3 in 1:(l_small + 1)^2
        gs = [gaunt_coefficients(ls[n2],ms[n2],ls[n3],ms[n3],ls[i],ms[i]) for i in eachindex(ls)]

        @test 4pi * (-1)^(ms[n3]+ms[n2]) * (1.0im)^(ls[n3]-ls[n2]) * Ys[n3] * conj(Ys[n2]) ≈
            sum( (1.0im).^(-ls) .* (-1.0).^ms .* conj.(Ys) .* gs) atol = 1e-12
    end

end
