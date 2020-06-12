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

    normal = [-1.0,0.0] # an outward normal to the surface
    material = Material(Halfspace(normal),specie)
    psource(θin) = PlaneSource(medium, [cos(θin),sin(θin)])

    k_effs = wavenumbers(ω, medium, material.species;
        tol=tol, basis_order = basis_order,
        num_wavenumbers = 15)
        # , mesh_points = 30, mesh_size = 2.0)

    function Rerror(θ)

        wave_effs = [
            WaveMode(ω, k_eff, psource(θ), material; tol=tol, basis_order = basis_order)
        for k_eff in k_effs]

        match_ws = MatchPlaneWaveMode(ω, psource(θ), material;
            basis_order = basis_order,
            tol = tol, wave_effs = wave_effs,
            max_size = 900,
        );

        Rm = reflection_coefficient(ω, match_ws, psource(θ), material)
        Rw = wienerhopf_reflection_coefficient(ω, psource(θ), material;
                tol=tol,
                basis_order = basis_order
        )
        return abs(Rm-Rw)
    end

    Rerror(0.0) # 3.81062386937501e-5
    Rerror(0.5) # 3.850176856470023e-5
    Rerror(1.4) # 0.00039393481919100916

    @test Rerror(0.0) < 1e-4
    @test Rerror(0.5) < 1e-4
    @test Rerror(1.4) < 1e-3

end

# Note the wiener-hopf method is inaccurate or takes along time to calculate a very low-frequency solution
@testset "compare wienger-hopf for low-ish frequency" begin

    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    specie = Specie(Acoustic(2; ρ=0.0, c=0.01), ms.Circle(0.01); volume_fraction = 0.3)
    s = ms.Circle(0.01)

    ω = 0.01
    basis_order=0

    tol = 1e-7

    eff_medium = effective_medium(medium, [specie])
    k_eff_low = ω / eff_medium.c

    k_eff_vol = wavenumber_low_volumefraction(ω, medium, [specie];
        basis_order = basis_order)

    normal = [-1.0,0.0]; # an outward normal to the surface
    material = Material(Halfspace(normal),specie);

    # Here we help the algorithm by giving some initial estimates. This is because k_eff_low does not work for either Dirichlet or when using basis_order = 0 (effective low frequency calculations assume basis_order = 1)

    k_effs = wavenumbers(ω, medium, material.species;
        tol=tol, basis_order = basis_order,
        box_k = [[-70.0,70.0],[0.0,50.0]] .* 6,
        num_wavenumbers = 5)

    function Rerror(θ)
        source = PlaneSource(medium, [cos(θ),sin(θ)]);

        wm = WaveMode(ω, k_effs[1], source, material; tol=tol, basis_order = basis_order, extinction_rescale = true)

        R1 = reflection_coefficient(ω, wm, source, material)
        R_low = reflection_coefficient(source, eff_medium, material.shape)
        Rw = wienerhopf_reflection_coefficient(ω, source, material;
                tol=tol,
                basis_order = basis_order,
                num_coefs = 20000 + 200*Int(round(100*θ))
        )
        return abs(2*R_low-Rw-R1)/abs(2*R_low)
    end

    @test Rerror(0.5) < 5e-4
    @test Rerror(1.0) < 5e-4
end
