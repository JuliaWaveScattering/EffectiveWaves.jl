using EffectiveWaves, Test

# this is for high volume and frequency
@testset "compare wienger-hopf and match method reflection" begin

    medium = Medium(1.0,1.0+0.0im)

    ω = 1.3

    specie = Specie(ρ=0.0, r=0.5, c=0.0, volfrac = 0.2)
    radius_multiplier = 1.01

    tol = 1e-9
    hankel_order=0
    t_vec = t_vectors(ω, medium, [specie]; hankel_order = hankel_order)

    ws = effective_waves(ω, medium, [specie];
        radius_multiplier = radius_multiplier,
        tol=tol, hankel_order = hankel_order,
        apply_meshing = true,
        num_wavenumbers = 20,
        mesh_points = 30, mesh_size = 2.);

    function Rerror(θ)

        match_ws = MatchWave(ω, medium, specie;
            radius_multiplier = radius_multiplier,
            hankel_order = hankel_order,
            tol = tol, wave_effs = ws[1:15],
            θin=θ,
            max_size = 600,
        );

        Rm = reflection_coefficient(ω, match_ws, medium, specie; θin = θ)
        Rw = wienerhopf_reflection_coefficient(ω, medium, [specie];
                θin = θ,
                tol=tol,
                radius_multiplier = radius_multiplier,
                hankel_order = hankel_order
        )
        return abs(Rm-Rw)/abs(Rw)
    end

    @test Rerror(0.0) < 1e-4
    @test Rerror(0.5) < 1e-4
    @test Rerror(1.4) < 1e-3
end

# Note the wiener-hopf method is inaccurate or takes along time to calculate a very low-frequency solution
@testset "compare wienger-hopf for low-ish frequency" begin

    medium = Medium(1.0,1.0+0.0im)
    ω = 0.01
    hankel_order=0

    # Always Dirichlet boundary conditions for hankel_order=0!
    specie = Specie(ρ = 0.0, r=0.01, c=6.0, volfrac = 0.3)

    tol = 1e-7
    t_vec = t_vectors(ω, medium, [specie]; hankel_order = hankel_order)

    eff_medium = effective_medium(medium, [specie])

    function Rerror(θ)
        ws = effective_waves(ω, medium, [specie];
            tol=tol, θin = θ,
            hankel_order = hankel_order,
            num_wavenumbers = 2);
        R1 = reflection_coefficient(ω, ws[1], medium, [specie]; θin = θ)
        R_low = reflection_coefficient_halfspace(medium, eff_medium; θin = θ)
        Rw = wienerhopf_reflection_coefficient(ω, medium, [specie];
                θin = θ,
                tol=tol,
                hankel_order = hankel_order,
                num_coefs = 20000 + 200*Int(round(100*θ))
        )
        return abs(2*R_low-Rw-R1)/abs(2*R_low)
    end

    @test Rerror(0.5) < 1e-3
    @test Rerror(1.0) < 1e-3
end
