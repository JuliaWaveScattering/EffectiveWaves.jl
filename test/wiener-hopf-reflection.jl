using EffectiveWaves, Test

# this is for high volume and frequency
@testset "compare wienger-hopf and match method reflection" begin

    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    specie = Specie(Particle(
        Acoustic(2; ρ=0.0, c=0.0), ms.Circle(0.5));
        volume_fraction = 0.2, exclusion_distance = 1.01
    )

    ω = 1.3

    tol = 1e-9
    basis_order=0

    ws = effective_wavemodes(ω, medium, [specie];
        tol=tol, basis_order = basis_order,
        apply_meshing = false,
        num_wavenumbers = 20,
        mesh_points = 30, mesh_size = 2.);

    function Rerror(θ)

        match_ws = MatchPlaneWaveMode(ω, medium, specie;
            basis_order = basis_order,
            tol = tol, wave_effs = ws[1:15],
            θin=θ,
            max_size = 600,
        );

        Rm = reflection_coefficient(ω, match_ws, medium, specie; θin = θ)
        Rw = wienerhopf_reflection_coefficient(ω, medium, [specie];
                θin = θ,
                tol=tol,
                basis_order = basis_order
        )
        return abs(Rm-Rw)/abs(Rw)
    end

    Rerror(0.0) # 3.771528229242978e-5
    Rerror(0.5) # 5.2531420000213484e-5
    Rerror(1.4) # 0.0005018041373008065

    @test Rerror(0.0) < 1e-4
    @test Rerror(0.5) < 1e-4
    @test Rerror(1.4) < 1e-3

end

# Note the wiener-hopf method is inaccurate or takes along time to calculate a very low-frequency solution
@testset "compare wienger-hopf for low-ish frequency" begin

    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering
    
    # Always Dirichlet boundary conditions for basis_order=0!
    specie = Specie(Particle(
        Acoustic(2; ρ=0.0, c=6.0), ms.Circle(0.01));
        volume_fraction = 0.3
    )

    ω = 0.01
    basis_order=0

    tol = 1e-7

    eff_medium = effective_medium(medium, [specie])

    function Rerror(θ)
        ws = effective_wavemodes(ω, medium, [specie];
            tol=tol, θin = θ,
            basis_order = basis_order,
            num_wavenumbers = 2,
            extinction_rescale = true);
        R1 = reflection_coefficient(ω, ws[1], medium, [specie]; θin = θ)
        R_low = reflection_coefficient(medium, eff_medium; θin = θ)
        Rw = wienerhopf_reflection_coefficient(ω, medium, [specie];
                θin = θ,
                tol=tol,
                basis_order = basis_order,
                num_coefs = 20000 + 200*Int(round(100*θ))
        )
        return abs(2*R_low-Rw-R1)/abs(2*R_low)
    end

    @test Rerror(0.5) < 1e-3
    @test Rerror(1.0) < 1e-3
end
