# Test that path and mesh methods find the same wavenumbers
using EffectiveWaves, Test

@testset "Mesh and path methods for wavenumbers" begin

    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    species = [
        Specie(Particle(Acoustic(2; ρ=8.0, c=1.1), ms.Circle(1.0)); volume_fraction=0.25),
    ]

    ω = 0.6;
    tol = 1e-6
    num_wavenumbers = 8;
    basis_order = 2

    k_effs_path = wavenumbers_path(ω, medium, species;
        mesh_points = 10,
        num_wavenumbers=num_wavenumbers,
        basis_order=basis_order,
        tol = tol
    )
    num_wavenumbers = min(num_wavenumbers, length(k_effs_path))

    k_effs_mesh = wavenumbers_mesh(ω, k_effs_path[1:num_wavenumbers], medium, species;
        basis_order=basis_order,
        tol = tol
    );

   @test maximum(abs.(k_effs_mesh[1:num_wavenumbers] - k_effs_path[1:num_wavenumbers])) < 10*tol

end
