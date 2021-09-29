

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
        exclusion_distance = separation_ratio
    );
    species = [s1]
    material = Material(Sphere(R),species);

    a12 = 2.0 * minimum(outer_radius(s) * s.exclusion_distance for s in material.species)

    tol = 1e-2

    Uinter = outgoing_translation_matrix(ω, medium, material;  basis_order = basis_order, tol = tol)

    rs = LinRange(a12,2R-a12,13);
    φs = LinRange(-pi,pi,13);
    θs = LinRange(0.0,pi,13);

    r2x = radial_to_cartesian_coordinates

    error = [
        begin
            U1 = outgoing_translation_matrix(medium, basis_order, ω, r2x([r,θ,φ]))
            norm(Uinter(r2x([r,θ,φ])) - U1) / norm(U1)
        end
    for r in rs, θ in θs, φ in φs];

    @test maximum(error) < 5 * tol

    @test mean(error) < tol

end
