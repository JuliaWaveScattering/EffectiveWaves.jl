@testset "Particle types and constructors" begin
    ω=0.4
    medium = Medium(ρ=10.,c=2.0+3.0im)
    p_dirichlet = Specie(0.0,0.1) # ρ=0.0, r =0.1
    p_neumann = Specie(Inf,0.1) # ρ=0.0, r =0.1
    Z_dirichlet = Zn(ω, p_dirichlet, medium, 0)
    Z_neumann   = Zn(ω, p_neumann, medium, 0)
    @test Z_dirichlet == besselj(0, 0.1*ω/medium.c)/hankelh1(0, 0.1*ω/medium.c)
    @test Z_neumann == besselj(1, 0.1*ω/medium.c)/hankelh1(1, 0.1*ω/medium.c)

    medium = Medium(ρ=0.,c=2.0+3.0im)
    @test_throws ErrorException Zn(ω, p_dirichlet, medium, 0)
    @test Z_neumann == Zn(ω, p_neumann, medium, 0)

    f_tmp(x,y; kws...) = 0.0
    @test gray_square(rand(10),rand(10), f_tmp) == 0.0
end
