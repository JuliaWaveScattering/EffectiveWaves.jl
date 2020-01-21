using EffectiveWaves, Test

@testset "match wave low volfrac" begin

using LinearAlgebra
## Low volume fraction
    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    exclusion_distance = 1.005

    specie = Specie(Particle(
        Acoustic(2; ρ=0.5, c=0.5), ms.Circle(0.4));
        volume_fraction=0.001,
        exclusion_distance=exclusion_distance
    )
    species = [specie]

    θin = 0.2
    tol = 1e-8
    basis_order = 2

    ωs = [0.2,1.2]

    wave_effs_arr = [
        effective_waves(ω, medium, [specie];
            basis_order=basis_order, tol = tol,  θin = θin,
            num_wavenumbers = 1
            #, mesh_points = 10, mesh_size = 2.0 #, max_Rek = 20.0, max_Imk = 20.0
            , extinction_rescale = false)
    for ω in ωs]

    k_eff_φs = wavenumber_low_volumefraction(ωs, medium, species; tol=tol,
        basis_order=basis_order)

    wave_eff_φs = [
        EffectiveWave(ωs[i], k_eff_φs[i], medium, species; tol = tol,
            basis_order=basis_order, θin = θin,
            extinction_rescale = true)
    for i in eachindex(ωs)]

    @test maximum(abs(wave_effs_arr[i][1].k_eff/k_eff_φs[i] - 1.0) for i in eachindex(ωs)) < 1e-7

    ds = [
        abs(dot(wave_eff_φs[i].amplitudes[:], wave_effs_arr[i][1].amplitudes[:]))
    for i in eachindex(ωs)]
    ds2 = [
        norm(wave_eff_φs[i].amplitudes[:])*norm(wave_effs_arr[i][1].amplitudes[:])
    for i in eachindex(ωs)]
    @test maximum(ds./ds2 .- 1.0) < tol

    match_ws = [
        MatchWave(ωs[i], medium, specie; θin = θin,
            max_size=150,
            tol = tol, wave_effs = wave_effs_arr[i][1:1])
    for i in eachindex(ωs)]

    @test maximum(match_error.(match_ws)) < 1e-6

    avg_wave_φs = [
        AverageWave(match_ws[i].discrete_wave.x, wave_eff_φs[i])
    for i in eachindex(ωs)]

    Rs_near = [reflection_coefficient(ωs[i], match_ws[i].discrete_wave, medium, specie) for i in eachindex(ωs)]
    Rs_near_φ = [reflection_coefficient(ωs[i], avg_wave_φs[i], medium, specie) for i in eachindex(ωs)]
    @test maximum(abs.(Rs_near_φ - Rs_near)) < 5e-6

    Rs = [reflection_coefficient(ωs[i], match_ws[i], medium, specie) for i in eachindex(ωs)]
    Rs_φ = [reflection_coefficient(ωs[i], wave_eff_φs[i], medium, species) for i in eachindex(ωs)]
    @test maximum(abs.(Rs_near_φ - Rs_near)) < 5e-6

    @test maximum(
        norm(match_ws[i].discrete_wave.amplitudes[:] .- avg_wave_φs[i].amplitudes[:])/norm(avg_wave_φs[i].amplitudes[:])
    for i in eachindex(ωs)) < 1e-3
end
