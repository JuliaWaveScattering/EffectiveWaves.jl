using EffectiveWaves, Test

@testset "match wave low volfrac" begin

## Low volume fraction
    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    exclusion_distance = 1.005

    specie = Specie(Particle(
        Acoustic(2; ρ=0.5, c=0.5), ms.Circle(0.4));
        volume_fraction=0.001,
        exclusion_distance=exclusion_distance
    )

    θin = 0.2
    tol = 1e-8
    basis_order = 2

    normal = [-1.0,0.0] # an outward normal to the surface
    material = Material(Halfspace(normal),[specie])
    source = PlaneSource(medium, [cos(θin),sin(θin)])

    ωs = [0.2,1.2]

    wave_effs_arr = map(ωs) do ω
        WaveModes(ω, source, material;
            basis_order=basis_order, tol = tol,
            extinction_rescale = false,
            # , mesh_points = 10, mesh_size = 2.0 #, max_Rek = 20.0, max_Imk = 20.0
            num_wavenumbers = 1)
    end

    k_eff_φs = wavenumber_low_volumefraction(ωs, medium, [specie]; basis_order=basis_order)

    wave_eff_φs = [
        WaveMode(ωs[i], k_eff_φs[i], source, material; tol = tol,
            basis_order=basis_order,
            extinction_rescale = true)
    for i in eachindex(ωs)]

    @test maximum(abs.([w[1].wavenumber for w in wave_effs_arr] - k_eff_φs)) < 1e-7

    # Check that the eigenmodes point the same direction
    ds = [
        abs(dot(wave_eff_φs[i].eigenvectors[:], wave_effs_arr[i][1].eigenvectors[:]))
    for i in eachindex(ωs)]
    ds2 = [
        norm(wave_eff_φs[i].eigenvectors[:])*norm(wave_effs_arr[i][1].eigenvectors[:])
    for i in eachindex(ωs)]
    @test maximum(ds./ds2 .- 1.0) < tol

    match_ws = [
        MatchPlaneWaveMode(ωs[i], source, material;
            max_size=150,
            tol = tol, wave_effs = wave_effs_arr[i][1:1])
    for i in eachindex(ωs)]

    @test maximum(match_error(m,material.shape) for m in match_ws) < 1e-6

    avg_wave_φs = [
        DiscretePlaneWaveMode(match_ws[i].discrete_wave.x, wave_eff_φs[i], material.shape)
    for i in eachindex(ωs)]

    Rs_near = [reflection_coefficient(ωs[i], match_ws[i].discrete_wave, source, material) for i in eachindex(ωs)]
    Rs_near_φ = [reflection_coefficient(ωs[i], avg_wave_φs[i], source, material) for i in eachindex(ωs)]
    @test maximum(abs.(Rs_near_φ - Rs_near)) < 5e-6

    Rs = [reflection_coefficient(ωs[i], match_ws[i], source, material) for i in eachindex(ωs)]
    Rs_φ = [reflection_coefficient(ωs[i], wave_eff_φs[i], source, material) for i in eachindex(ωs)]
    @test maximum(abs.(Rs_near_φ - Rs_near)) < 5e-6

    @test maximum(
        norm(match_ws[i].discrete_wave.amplitudes[:] .- avg_wave_φs[i].amplitudes[:])/norm(avg_wave_φs[i].amplitudes[:])
    for i in eachindex(ωs)) < 1e-3
end
