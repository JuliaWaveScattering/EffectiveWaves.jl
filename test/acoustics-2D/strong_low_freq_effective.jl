using EffectiveWaves, Test

# This case is numerically challenging, because wavenumber() has many roots close together. Make sure spacing in ω is small to help the optimisation method
@testset "strong scatterers and low frequency" begin
    medium = Acoustic(2; ρ=1.0, c=1.0)
    spatial_dim = 2

    ms = MultipleScattering # just in case Circle conflicts with a definition from another package.

    s1 = Specie(
        Acoustic(spatial_dim; ρ=5.0, c=0.002),ms.Circle(0.004);
        volume_fraction=0.2
    )
    s2 = Specie(
        Acoustic(spatial_dim; ρ=0.3, c=0.01),ms.Circle(0.002);
        volume_fraction=0.1
    )
    species = [s1, s2]
    tol=1e-7
    # ωs = [0.001,0.003]
    ωs = [0.001]
    eff_medium = effective_medium(medium, species)
    k_eff_lows = ωs./eff_medium.c

    k_eff_φs = wavenumber_low_volumefraction(ωs, medium, species; basis_order=1)

    k_effs_arr = [
        wavenumbers(ω, medium, species; tol=tol, num_wavenumbers=1, basis_order=1)
        # num_wavenumbers = 1 usually finds the wavenubmer with the smallest attenuation
    for ω in ωs]

    inds = [argmin(abs.(k_effs_arr[i] .- k_eff_φs[i])) for i in eachindex(ωs)]
    k_effs2 = [k_effs_arr[i][inds[i]] for i in eachindex(inds)]

    @test norm(k_effs2 - k_eff_lows)/norm(k_effs2) < 5e-7
    @test norm(k_effs2[1] - k_eff_lows[1])/norm(k_eff_lows[1]) < 5e-7
    @test norm(k_effs2 - k_eff_φs)/norm(k_effs2) < 0.01
end
