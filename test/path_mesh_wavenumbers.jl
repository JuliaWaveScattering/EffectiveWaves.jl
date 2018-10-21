# Test that path and mesh methods find the same wavenumbers
using EffectiveWaves, Test

@testset "Mesh and path methods for wavenumbers" begin

    # medium = WaterDistilled
    medium = Medium(1.0,1.0+0.0im)

    species = [
        Specie(ρ=8.0, r=1.0, c=1.1, volfrac=0.25)
    ];

    ω = 0.6;
    tol = 1e-8
    num_wavenumbers = 8;
    hankel_order = 2

    k_effs_path = wavenumbers_path(ω, medium, species;
        mesh_size = 2.0,
        mesh_points = 10,
        num_wavenumbers=num_wavenumbers+4,
        hankel_order=hankel_order,
        tol = tol
    )

    k_effs_mesh = wavenumbers_mesh(ω, k_effs_path[1:num_wavenumbers], medium, species;
        hankel_order=hankel_order,
        tol = tol
    );

   @test maximum(abs.(k_effs_mesh[1:num_wavenumbers] - k_effs_path[1:num_wavenumbers])) < tol

end
