using EffectiveWaves, Test
using LinearAlgebra

@testset "3D pair-correlation" begin

    # NOTE: the package has only Percus-Yevick and hole-correction implemented
    r = 1.0

    # choose the type of pair correlation
    pair_type = PercusYevick(3; rtol = 1e-2, maxsize = 50)

    s = Specie(
        Acoustic(3; ρ = 10.0, c = 10.0),
        Sphere(r),
        volume_fraction = 0.3,
        exclusion_distance = 1.01
    );

    micro = Microstructure(s, pair_type);

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
