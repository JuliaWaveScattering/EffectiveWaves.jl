using EffectiveWaves, Test

# large volume fraction scatterers,  small size amd on strong scatterer. This is a problamatic case.
@testset "large volume fraction and low frequency" begin
    # should pass for the list of angular frequencies below
    # ωs = [0.001, 2.0, 9.0]
    ω = 0.001
    basis_order = 1

    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    s1 = Specie(Acoustic(2; ρ=5.0, c=1.2),ms.Circle(0.004); volume_fraction=0.4)
    s2 = Specie(Acoustic(2; ρ=0.3, c=0.4),ms.Circle(0.002); volume_fraction=0.3)

    species = [s1, s2]

    eff_medium = effective_medium(medium, species)

    tol = 1e-5
    k_eff_low = ω/eff_medium.c
    k_effs = wavenumbers(ω, medium, species; tol = tol, num_wavenumbers=1, basis_order=basis_order)
    i = findmin(abs.(k_effs .- k_eff_low))[2]
    k_eff = k_effs[i]

    @test abs(k_eff - k_eff_low)/norm(k_eff_low) < tol

    # let's now assume the particles are fill a halfspace. Then we can calculate a reflection coefficent from this halfspace
    normal = [-1.0,0.0] # an outward normal to the surface
    material = Material(Halfspace(normal),species)

    # define a plane wave source travelling at a 45 degree angle in relation to the material
    source = PlaneSource(medium, [cos(pi/4.0),sin(pi/4.0)])

    R = begin
        wave = WaveMode(ω, k_eff, source, material; basis_order=basis_order)
        reflection_coefficient(ω, wave, source, material)
    end
    R_low2 = begin
        wave2 = WaveMode(ω, k_eff_low, source, material; basis_order=basis_order)
        reflection_coefficient(ω, wave2, source, material)
    end
    R_low = reflection_transmission_coefficients(ω, source, eff_medium, material.shape)[1]

    @test norm(R_low - R) < tol
    @test norm(R_low2 - R) < tol
end
