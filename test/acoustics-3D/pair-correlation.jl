# NOTE: when checking against the literature, we find the same values for the PercusYevick approxiation, but do not find the same values for our method for the MonteCarloPairCorrelation

using EffectiveWaves, Test
using LinearAlgebra

@testset "3D pair-correlation" begin

    # NOTE: the package has only Percus-Yevick, MonteCarloPairCorrelation, and HoleCorrection implemented

    # choose the type of pair correlation
    # pairtype = PercusYevick(3; rtol = 1e-3, meshsize = 0.05, maxlength = 250)
    pairtype = PercusYevick(3; rtol = 1e-3, meshsize = 0.05, maxlength = 50)

    pairtype_mc = MonteCarloPairCorrelation(3; rtol = 1e-3, maxlength = 50, iterations = 1)

    # pmcs = [
    #     MonteCarloPairCorrelation(3; rtol = 1e-3, maxlength = 50, meshsize = 0.1, iterations = 20, numberofparticles = 60000),
    #     MonteCarloPairCorrelation(3; rtol = 1e-3, maxlength = 120, meshsize = 0.02,iterations = 200),
    #     MonteCarloPairCorrelation(3; rtol = 1e-3, maxlength = 60, iterations = 2000, numberofparticles = 2000)
    # ]

    s = Specie(
        Acoustic(3; ρ = 10.0, c = 10.0),
        # Sphere(1.5),
        Sphere(0.5),
        volume_fraction = 0.2,
        separation_ratio = 1.0
    );

    micro_mc = Microstructure(s, pairtype_mc);

    # Note that the volume fraction achieved with the MonteCarlo approach is not exactly the same as the one requested
    vol_mc = micro_mc.paircorrelations[1].number_density * volume(s)

    # PercusYevick only seems to match when using a lower volume fraction. We have yet to determine why.

    s_py = Specie(
        Acoustic(3; ρ = 10.0, c = 10.0),
        Sphere(0.5),
        volume_fraction = vol_mc,
        # volume_fraction = vol_mc*0.76,  # <--- more accurate, but don't know why
        separation_ratio = 1.0
    );
    micro = Microstructure(s_py, pairtype);

    # gs = [micro.paircorrelations[1], micro.paircorrelations[1]]

    ω = 1.2
    basis_order = 1
    basis_field_order = 3
    medium = Acoustic(3; ρ=1.0, c=1.0)

    kps_py = wavenumbers(ω, medium, micro;
        basis_order = basis_order, num_wavenumbers = 4
    )
    kps = wavenumbers(ω, medium, s_py;
        basis_order = basis_order, num_wavenumbers = 4
    )

    psource = PlaneSource(medium, [0.0,0.0,1.0]);

    # choose the size and position of the spherical domain of the material
    R = 10.0
    material = Material(medium,Sphere(R),micro);

    @test Symmetry(material,psource) == AzimuthalSymmetry{3}()

    # for the analytic solution, we need the wavemode first
    wave_py = WaveMode(ω, kps_py[1], psource, material;
        basis_order = basis_order,
        basis_field_order = basis_field_order
    )

    wave = WaveMode(ω, kps_py[1], psource, Material(medium,Sphere(R),[s_py]);
        basis_order = basis_order,
        basis_field_order = basis_field_order
    )

    # and then the outward spherical wave coefficients
    scat_coefficients_py = material_scattering_coefficients(wave_py)
    scat_coefficients = material_scattering_coefficients(wave)


## Test for a sphere and compare with the integral method

end


## NOTE to get accurate comparison between MC and PY, for the code here,

# using Plots
#
# # for vol = 5%
# # micro_mc5 = deepcopy(micro_mc)
# plot(micro.paircorrelations[1].r, micro.paircorrelations[1].dp)
# plot!(micro_mc.paircorrelations[1].r, micro_mc.paircorrelations[1].dp * 1.42,
# linestyle = :dash)
#
# # for vol = 10%
# plot(micro.paircorrelations[1].r, micro.paircorrelations[1].dp)
# plot!(micro_mc.paircorrelations[1].r, micro_mc.paircorrelations[1].dp * 1.38)
#
# # for vol = 15%
# plot(micro.paircorrelations[1].r, micro.paircorrelations[1].dp)
# plot!(micro_mc.paircorrelations[1].r, micro_mc.paircorrelations[1].dp * 1.38)
#
# # for vol = 20%
# plot(micro.paircorrelations[1].r, micro.paircorrelations[1].dp)
# plot!(micro_mc.paircorrelations[1].r, micro_mc.paircorrelations[1].dp * 1.38)
