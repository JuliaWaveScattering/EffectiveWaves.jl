# Test that path and mesh methods find the same wavenumbers
using EffectiveWaves, Test

@testset "bisection and path methods for wavenumbers" begin

    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    ω = 0.6;
    tol = 1e-6
    basis_order = 2
    
    species = [
        Specie(Particle(Acoustic(2; ρ=8.0, c=1.1), ms.Circle(1.0)); volume_fraction=0.25),
    ]
    micro = Microstructure(medium,species)

    num_wavenumbers = 1;
    k_effs_path = wavenumbers_path(ω, micro;
        mesh_points = 10,
        num_wavenumbers=num_wavenumbers,
        basis_order=basis_order,
        tol = tol
    )
    
    k_effs_bi = wavenumbers_bisection_robust(ω, micro;
        num_wavenumbers=num_wavenumbers,
        basis_order=basis_order,
        tol = tol
    )
    @test maximum(abs.(k_effs_bi[1] - k_effs_path[1])) < 10*tol

    ωs = LinRange(0.1, ω, 30)
    k_effs = wavenumbers(ωs, micro; tol = tol, basis_order = basis_order)
    @test abs(k_effs[end] - k_effs_bi[1]) < 10*tol

    num_wavenumbers = 8;

    k_effs_path = wavenumbers_path(ω, micro;
        mesh_points = 20,
        mesh_size = 2,
        num_wavenumbers=num_wavenumbers,
        basis_order=basis_order,
        tol = tol
    )
    
    k_effs_bi = wavenumbers_bisection_robust(ω, micro;
        num_wavenumbers=num_wavenumbers,
        basis_order=basis_order,
        tol = tol
    )

    @test maximum(abs.(k_effs_bi[1:5] - k_effs_path[1:5])) < 10*tol

## Void species[]
    species = [
        Specie(Particle(Acoustic(2; ρ=0.8, c=0.1), ms.Circle(1.0)); volume_fraction=0.25),
    ]
    micro = Microstructure(medium,species)

    num_wavenumbers = 1;
    k_effs_path = wavenumbers_path(ω, micro;
        num_wavenumbers=num_wavenumbers,
        basis_order=basis_order,
        tol = tol
    )
    
    k_effs_bi = wavenumbers_bisection_robust(ω, micro;
        num_wavenumbers=num_wavenumbers,
        basis_order=basis_order,
        tol = tol
    )
    @test maximum(abs.(k_effs_bi[1] - k_effs_path[1])) < 10*tol

    # Case below doesn't work. Because the wavenumber with the smaller imaginary part for ω is not on the same branch as the wavenumber with the smaller imaginary part for lower frequenices. This can be seen be running the below.
    # ωs = LinRange(0.6, 0.01, 5200)
    # k_effs = wavenumbers(ωs, micro; tol = tol, basis_order = basis_order)
    # @test abs(k_effs[end] - k_effs_bi[1]) < 10*tol

    num_wavenumbers = 8;

    k_effs_path = wavenumbers_path(ω, micro;
        mesh_points = 10,
        num_wavenumbers=num_wavenumbers,
        basis_order=basis_order,
        tol = tol
    )
    
    k_effs_bi = wavenumbers_bisection_robust(ω, micro;
        num_wavenumbers=num_wavenumbers,
        basis_order=basis_order,
        tol = tol
    )

    @test maximum(abs.(k_effs_bi[1:6] - k_effs_path[1:6])) < 10*tol

    # using Plots
    # scatter(real.(k_effs_path),imag.(k_effs_path))
    # # plot!(xlims = (0.4, 0.8), ylims = (-0.01, 0.3))
    
    # scatter!(real.(k_effs_bi),imag.(k_effs_bi))
    # num_wavenumbers = min(num_wavenumbers, length(k_effs_path))

end
