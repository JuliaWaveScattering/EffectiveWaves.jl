

@testset "Test interpolation of outgoing_translation_matrix" begin

    particle_medium = Acoustic(3; ρ=10.0, c=80.0);
    medium = Acoustic(3; ρ=1.0, c=1.0);

    R = 3.0
    r = 1.0

    separation_ratio = 1.02
    ka = 0.1
    k = ka / r
    vol_fraction = 0.3

    basis_order = 2

    ω = k * real(medium.c)

    s1 = Specie(
        particle_medium, Sphere(r);
        number_density = vol_fraction / volume(Sphere(r)),
        separation_ratio = separation_ratio
    );
    species = [s1]
    material = Material(Sphere(R),species);

    a12 = 2.0 * minimum(outer_radius(s) * s.separation_ratio for s in material.microstructure.species)

    tol = 1e-2

    Uinter = outgoing_translation_matrix(ω, material;  basis_order = basis_order, tol = tol)

    rs = LinRange(a12,2R-a12,13);
    φs = LinRange(-pi,pi,13);
    θs = LinRange(0.0,pi,13);

    r2x = radial_to_cartesian_coordinates

    error = [
        begin
            U1 = outgoing_translation_matrix(medium, basis_order, basis_order, ω, r2x([r,θ,φ]))
            norm(Uinter(r2x([r,θ,φ])) - U1) / norm(U1)
        end
    for r in rs, θ in θs, φ in φs];

    @test maximum(error) < 5 * tol

    @test mean(error) < tol

end

@testset "Legendre polynomial approximation" begin

    P = Legendre{Float64}()

    polynomial_order = 25
    mesh_points = 3*(polynomial_order + 1)
    # mesh_points = (polynomial_order + 1)^2

    rs = LinRange(0,5,mesh_points)
    r2s = LinRange(0,5,4*mesh_points)

    r1_max = maximum(rs);
    rbars = 2.0 .* rs ./ r1_max .- (1.0);
    r2bars = 2.0 .* r2s ./ r1_max .- (1.0);

    data = cos.(rs) .* (1.0 .+ 4 .* exp.(-10 .* rs.^4))
    data2 = cos.(r2s) .* (1.0 .+ 4 .* exp.(-10 .* r2s.^4))
    # plot(rs,data)
    # plot!(r2s,data2)

    Pmat = P[rbars, 1:(polynomial_order + 1)];
    P2mat = P[r2bars, 1:(polynomial_order + 1)];

    # pls_arr * transpose(Pmat) ~ data
    # Pmat * transpose(pls_arr) ~ transpose(data)
    pls_arr = Pmat \ data

    ws = (0:polynomial_order) .+ 1/2
    σs = integration_scheme(rbars; scheme = :trapezoidal)
    tls_arr = ws .* (transpose(Pmat) * (σs .* data))

    λ = 1e-6
    cls_arr = inv(transpose(Pmat) * Pmat + λ*I) * transpose(Pmat) * data

    # norm(Pmat * pls_arr - data) / norm(data)
    # norm(Pmat * tls_arr - data) / norm(data)
    # norm(Pmat * cls_arr - data) / norm(data)
    @test norm(P2mat * pls_arr - data2) / norm(data2) < 5e-3
    @test norm(P2mat * tls_arr - data2) / norm(data2) < 0.6
    @test norm(P2mat * cls_arr - data2) / norm(data2) < 5e-3

# using Plots
#
#     scatter(rs,data)
#     plot!(rs,Pmat * pls_arr, linestyle = :dash)
#     plot!(r2s,P2mat * pls_arr, linestyle = :dash)
#     plot!(r2s,P2mat * tls_arr, linestyle = :dash)
#     plot!(r2s,P2mat * cls_arr, linestyle = :dash)

end
