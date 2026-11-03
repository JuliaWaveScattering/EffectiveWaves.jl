using EffectiveWaves, Test, LinearAlgebra

@testset "3D planar symmetry with particles much smaller than the wavelength" begin

    # A halfspace holding both "fines", with k a between 6e-5 and 6e-4, and "coarse"
    # particles, with k a between 0.004 and 0.08. The columns of the eigensystem carry the
    # number density n ~ a⁻³ of their specie and the rows its T-matrix T ~ (k a)³, so the
    # fines' components of the eigenvector are O(T) relative to the coarse ones. An SVD of
    # the eigensystem as it stands cannot resolve them, and it also returns spurious
    # near-null vectors made of the fines alone, which made `WaveMode` find "more than one
    # eigenvector" and fail with a BoundsError. `eigenvectors` now equilibrates the rows
    # and columns of the eigensystem before its SVD when their norms are this spread out.

    medium = Acoustic(3; ρ = 1.0, c = 1.0)
    particle = Acoustic(3; ρ = 3.96, c = 6.67)
    psource = PlaneSource(medium, [0.0, 0.0, 1.0])
    halfspace = Halfspace([0.0, 0.0, -1.0])

    fines = [
        Specie(particle, Sphere(a); volume_fraction = 0.001)
    for a in exp.(range(log(0.003), log(0.03), length = 10))]
    coarse = [
        Specie(particle, Sphere(a); volume_fraction = 0.002)
    for a in exp.(range(log(0.2), log(4.0), length = 10))]
    species = [fines; coarse]
    material = Material(medium, halfspace, species)

    ω = 0.02
    basis_order = 3
    k_eff = wavenumber_low_volumefraction(ω, medium, species; basis_order = basis_order)

    # A single eigenvector, both at this wavenumber and at one 0.01% away from it, as a
    # wavenumber from an approximation would be
    for k in (k_eff, k_eff * (1 + 1e-4))
        F = eigenvectors(ω, k, psource, material; basis_order = basis_order, tol = 1e-2)
        @test size(F, 3) == 1
    end

    wavemode = WaveMode(ω, k_eff, psource, material; basis_order = basis_order)
    R = reflection_coefficient(wavemode, psource, material)

    # At this low frequency |R| is that of the long-wavelength effective medium. Its phase
    # is not, as the particle centres are excluded from a layer next to the boundary.
    R_low = reflection_transmission_coefficients(ω, psource, effective_medium(medium, species), halfspace)[1]
    @test abs(R) ≈ abs(R_low) rtol = 1e-2
end
