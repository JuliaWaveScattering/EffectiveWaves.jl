# Test that path and mesh methods find the same wavenumbers
using EffectiveWaves, Test

@testset "bisection and path methods for wavenumbers" begin

    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    species = [
        Specie(Particle(Acoustic(2; ρ=8.0, c=1.1), ms.Circle(1.0)); volume_fraction=0.25),
    ]

    ω = 0.6;
    tol = 1e-6
    num_wavenumbers = 8;
    basis_order = 2

    micro = Microstructure(medium,species)
    k_effs_path = wavenumbers_path(ω, micro;
        mesh_points = 10,
        num_wavenumbers=num_wavenumbers,
        basis_order=basis_order,
        tol = tol
    )
    
    using Plots
    scatter(real.(k_effs_path),imag.(k_effs_path))
    plot!(xlims = (0.4, 0.8), ylims = (-0.01, 0.3))
    # num_wavenumbers = min(num_wavenumbers, length(k_effs_path))

    # Create grids for real and imaginary parts of kz
    kz_real_vals = -0.4:0.01:1.0
    kz_imag_vals = -0.4:0.01:0.3

    kz_real_vals = -5.4:0.01:5.0
    kz_imag_vals = -0.4:0.01:3.3

     # Calculate disp values for all combinations
    Z_real = zeros(length(kz_imag_vals), length(kz_real_vals))
    Z_imag = zeros(length(kz_imag_vals), length(kz_real_vals))
    
    disp = dispersion_complex(ω, micro)

    for (i, kz_imag) in enumerate(kz_imag_vals)
        for (j, kz_real) in enumerate(kz_real_vals)
            kz = kz_real + im*kz_imag
            disp_val = disp(kz)
            Z_real[i, j] = real(disp_val)
            Z_imag[i, j] = imag(disp_val)
        end
    end
    
    # Plot heatmaps
    using Plots
    p1 = heatmap(kz_real_vals, kz_imag_vals, Z_real,
                 xlabel="Re(kz)", ylabel="Im(kz)",
                 title="Real part of disp(ω=$ω, kz)",
                 c=:balance, clims = (-2,2))
    
    p2 = heatmap(kz_real_vals, kz_imag_vals, Z_imag,
                 xlabel="Re(kz)", ylabel="Im(kz)",
                 title="Imaginary part of disp(ω=$ω, kz)",
                 c=:balance, clims = (-1,1))
    
    plot(p1, p2, layout=(2,1), size=(1200, 500))

    # k_effs_mesh = wavenumbers_mesh(ω, k_effs_path[1:num_wavenumbers], micro;
    #     basis_order=basis_order,
    #     tol = tol
    # );

   @test maximum(abs.(k_effs_mesh[1:num_wavenumbers] - k_effs_path[1:num_wavenumbers])) < 10*tol

end
