using Test, LinearAlgebra

@testset "Tests for complex transmission wavevector" begin

    k_eff = rand(-1:0.1:1.0) + rand(-1:0.1:1.0) * im
    incident_wavevector = rand(-1:0.1:1.0,3) + rand(-1:0.1:1.0,3) .* im
    surface_normal = rand(-1:0.1:1.0,3)
    surface_normal = surface_normal / norm(surface_normal)

    k_vec = transmission_wavevector(k_eff, incident_wavevector, surface_normal)

    @test sum(x^2 for x in k_vec) ≈ k_eff^2
    @test imag(dot(-surface_normal, k_vec)) > 0

    k_eff = rand(-1:0.1:1.0) + rand(-1:0.1:1.0) * im
    incident_wavevector = rand(-1:0.1:1.0,2) + rand(-1:0.1:1.0,2) .* im
    surface_normal = rand(-1:0.1:1.0,2)
    surface_normal = surface_normal / norm(surface_normal)

    k_vec = transmission_wavevector(k_eff, incident_wavevector, surface_normal)

    @test sum(x^2 for x in k_vec) ≈ k_eff^2
    @test imag(dot(-surface_normal, k_vec)) > 0

    k_eff = rand(-1:0.1:1.0) + rand() * im
    k = rand()
    θin = rand(-pi/2.0:0.2:pi/2.0)
    θeff = transmission_angle(k+0.0im, k_eff, θin)

    k_vec2 = k_eff .* [cos(θeff+pi/2.0),sin(θeff+pi/2.0)]

    surface_normal = [0.0,-1.0];
    incident_wavevector = k*[cos(θin+pi/2.0),sin(θin+pi/2.0)]
    k_vec = transmission_wavevector(k_eff, incident_wavevector, surface_normal)

    @test norm(k_vec - k_vec2) < 1e-10

end
