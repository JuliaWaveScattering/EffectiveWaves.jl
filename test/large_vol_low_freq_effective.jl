using EffectiveWaves, Test

# large volume fraction scatterers,  small size amd on strong scatterer. This is a problamatic case.
@testset "large volume fraction and low frequency" begin
    # should pass for the list of angular frequencies below
    # ωs = [0.001, 2.0, 9.0]
    ω = 0.001

    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    p1 = Particle(Acoustic(2; ρ=5.0, c=1.2),ms.Circle(0.004))
    p2 = Particle(Acoustic(2; ρ=0.3, c=0.4),ms.Circle(0.002))

    species = [
        Specie(p1; volume_fraction=0.4),
        Specie(p2; volume_fraction=0.3)
    ]

    eff_medium = effective_medium(medium, species)

    tol = 1e-6
    k_eff_low = ω/eff_medium.c
    k_effs = wavenumbers(ω, medium, species; tol = tol, num_wavenumbers=1, basis_order=1)
    i = findmin(abs.(k_effs .- k_eff_low))[2]
    k_eff = k_effs[i]

    @test abs(k_eff - k_eff_low)/norm(k_eff_low) < 10*tol

    # let's now assume the particles are fill a halfspace. Then we can calculate a reflection coefficent from this halfspace
    normal = [-1.0,0.0] # an outward normal to the surface
    material = Material(Halfspace(normal),species)

    # define a plane wave source travelling directly towards the material
    source = PlaneSource(medium, -normal)

    R = begin
        wave = effective_wavemode(ω, k_eff, source, material)
        reflection_coefficient(ω, wave, source, material)
    end
    R_low2 = begin
        wave = effective_wavemode(ω, k_eff_low, source, material)
        reflection_coefficient(ω, wave, source, material)
    end
    R_low = reflection_coefficient(medium, eff_medium)

    @test norm(R_low - R) < tol
    @test norm(R_low2 - R) < tol
end
