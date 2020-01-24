# This package deals with complex spherical coordinates and complex spherical harmonics. Many of the identities used are not well known, so we pay special attention to tests for these complex functions.

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


@testset "Tests for complex transmission angle" begin

    surface_normal = [0.0,0.0,-1.0]
    xs = [rand(-1.01:0.1:1.0,3) + rand(-1.01:0.1:1.0,3)*im for i = 1:100]
    @test maximum(abs(cartesian_to_radial_coordiantes(x)[2] - transmission_angle(x, surface_normal)) for x in xs) < 1e-14

    surface_normal = [-1.0,0.0]
    xs = [rand(-1.01:0.1:1.0,2) + rand(-1.01:0.1:1.0,2)*im for i = 1:100]
    @test maximum(abs(cartesian_to_radial_coordiantes(x)[2] - transmission_angle(x, surface_normal)) for x in xs) < 1e-14

    N = 10
    is = 1:N
    ns = [rand(-1.01:0.1:1.0,2) for i = is]
    ns = ns ./ norm.(ns)

    v_ins = [rand(-1.01:0.1:1.0,2) for i = is]
    θ_ins = transmission_angle.(v_ins, ns)

    k_effs = rand(-1.01:0.1:1.0,N) + rand(-1.0:0.1:1.0,N) .* im
    v_effs = transmission_wavevector.(k_effs, v_ins, ns)
    θ_effs = transmission_angle.(v_effs, ns)

    # Check snells law, i.e. that the components orthogonal to the surface are the same for the two wavevectors
    @test maximum(norm(v_effs[i] - v_ins[i] - ns[i] .* dot(ns[i], v_effs[i] - v_ins[i])) for i in is) < 1e-14
    @test maximum(abs.((k_effs  .* sin.(θ_effs)) .^2 - (norm.(v_ins)  .* sin.(θ_ins)) .^2)) < 1e-14

    # Check that the effective wave attenuates as it moves into the material
    @test all(imag.(dot.(- ns,v_effs)) .> 0)

    # NOTE sum(v_effs[i] .^2) ≈ k_effs[i] ^2, however sqrt(sum(v_effs[i] .^2)) is not necessarily equal to k_effs[i] as they may be on different branches

    # When chooseing k_effs that have an imaginary part that makes them decay when transmitting into the surface, then we can recover k_effs from the v_effs
    ns = [ -rand(0.2:0.1:1.0,2) for i = is]
    ns = ns ./ norm.(ns)

    v_ins = [rand(0.2:0.1:1.0,2) for i = is]

    θ_ins = transmission_angle.(v_ins, ns)

    # choose positive imaginary part
    k_effs = rand(-1.01:0.1:1.0,N) + rand(0.1:0.1:1.0,N) .* im
    v_effs = transmission_wavevector.(k_effs, v_ins, ns)
    θ_effs = transmission_angle.(v_effs, ns)

    @test maximum(abs.(k_effs  .* sin.(θ_effs) - norm.(v_ins)  .* sin.(θ_ins)) ) < 1e-14

    axises = [ [-1.0,0.0] for i = is]
    θ_ns = transmission_angle.(-ns, axises)

    [norm(v_ins[i] - norm(v_ins[i]) .* [cos(θ_ins[i] + θ_ns[i]), sin(θ_ins[i] + θ_ns[i])]) for i in is]
    [norm(v_effs[i] - k_effs[i] .* [cos(θ_effs[i] + θ_ns[i]), sin(θ_effs[i] + θ_ns[i])]) for i in is]

 k_effs .* sin.(θ_effs .+ pi/2.0) - [v_effs[i][2] for i in is]

# test for 3-dimensions

    ns = [rand(-1.01:0.1:1.0,3) for i = is]
    ns = ns[is] ./ norm.(ns[is])

    v_ins = [rand(-1.01:0.1:1.0,3) for i = is]
    θ_ins = transmission_angle.(v_ins, ns)

    v_effs = transmission_wavevector.(k_effs, v_ins, ns)
    θ_effs = transmission_angle.(v_effs, ns)

    # Check snells law, i.e. that the components orthogonal to the surface are the same for the two wavevectors
    @test maximum(norm(v_effs[i] - v_ins[i] - ns[i] .* dot(ns[i], v_effs[i] - v_ins[i])) for i in is) < 1e-14
    @test maximum(abs.((k_effs  .* sin.(θ_effs)) .^2 - (norm.(v_ins)  .* sin.(θ_ins)) .^2)) < 1e-14

    # Check that the effective wave attenuates as it moves into the material
    @test all(imag.(dot.(- ns,v_effs)) .> 0)
end
