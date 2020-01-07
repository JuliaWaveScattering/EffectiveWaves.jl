# Test that path and mesh methods find the same wavenumbers
using EffectiveWaves, Test

@testset "Mesh and path methods for wavenumbers" begin

    medium = Medium(ρ=1.0, c=1.0)

    species = [
        Specie(ρ=8.0, r=1.0, c=1.1, volfrac=0.25)
    ];

    ω = 0.6;
    tol = 1e-9
    num_wavenumbers = 8;
    hankel_order = 2
    dim = 2

    k_effs_path = wavenumbers_path(ω, medium, species;
        mesh_size = 2.0,
        mesh_points = 10,
        max_Imk = 4.0,
        dim = dim,
        num_wavenumbers=num_wavenumbers,
        hankel_order=hankel_order,
        tol = tol
    )

    k_effs_mesh = wavenumbers_mesh(ω, k_effs_path[1:num_wavenumbers], medium, species;
        dim = dim,
        hankel_order=hankel_order,
        tol = tol
    );

   @test maximum(abs.(k_effs_mesh[1:num_wavenumbers] - k_effs_path[1:num_wavenumbers])) < 10*tol

end
