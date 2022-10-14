using EffectiveWaves, Test
using LinearAlgebra

@testset "3D pair-correlation" begin

    # NOTE: the package has only Percus-Yevick, MonteCarloPairCorrelation, and HoleCorrection implemented

    # choose the type of pair correlation
    pairtype = PercusYevick(3; rtol = 1e-3, meshsize = 0.05, maxlength = 250)

    pairtype = PercusYevick(3; rtol = 1e-3, meshsize = 0.05, maxlength = 250)
    pairtype_mc = MonteCarloPairCorrelation(3; rtol = 1e-3, meshsize = 0.09, maxlength = 250, iterations = 150)
    # pairtype = MonteCarloPairCorrelation(3; rtol = 1e-3, meshsize = 0.09, maxlength = 250, iterations = 150)

    s = Specie(
        Acoustic(3; ρ = 10.0, c = 10.0),
        Sphere(1.0),
        volume_fraction = 0.15,
        seperation_ratio = 1.0
    );

    # Using Monte-carlo is far heavier
    micro_mc = Microstructure(s, pairtype_mc);

    # Note that the volume fraction achieved with the MonteCarlo approach is not exactly the same as the one requested
    micro_mc.paircorrelations[1].number_density * volume(s)

    # s = Specie(
    #     Acoustic(3; ρ = 10.0, c = 10.0),
    #     Sphere(1.0),
    #     volume_fraction = 0.123 ,
    #     seperation_ratio = 1.0
    # );
    micro = Microstructure(s, pairtype);

    using Plots

    plot(micro.paircorrelations[1].r, micro.paircorrelations[1].dp)
    plot!(micro_mc.paircorrelations[1].r, (micro_mc.paircorrelations[1].dp .+ 1) .* 0.96 .- 1.0)

    length(micro.paircorrelations[1].dp)

    ω = 0.4
    basis_order = 1
    basis_field_order = 3
    medium = Acoustic(3; ρ=1.0, c=1.0)

    kps = wavenumbers(ω, medium, micro;
        basis_order = basis_order, num_wavenumbers = 4
    )

    psource = PlaneSource(medium, [0.0,0.0,1.0]);

    # choose the size and position of the spherical domain of the material
    R = 10.0
    material = Material(Sphere(R),micro);

    @test Symmetry(material,psource) == AzimuthalSymmetry{3}()

    # for the analytic solution, we need the wavemode first
    wave = WaveMode(ω, kps[1], psource, material;
        basis_order = basis_order,
        basis_field_order = basis_field_order
    )

    # and then the outward spherical wave coefficients
    scat_coefficients = material_scattering_coefficients(wave)


end
