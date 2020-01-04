using EffectiveWaves, Test

@testset "match wave low volfrac" begin

using LinearAlgebra
## Low volume fraction
    medium = Medium(1.0,1.0+0.0im)
    specie = Specie(ρ=0.5, r=0.4, c=0.5, volfrac=0.001)
    species = [specie]
    θin = 0.2
    tol = 1e-9
    hankel_order = 2

    ωs = [0.2,1.2]
    radius_multiplier = 1.005

    wave_effs_arr = [
        effective_waves(ω, medium, [specie];
            hankel_order=hankel_order, tol = tol,  θin = θin,
            radius_multiplier = radius_multiplier,
            num_wavenumbers = 1
            #, mesh_points = 10, mesh_size = 2.0 #, max_Rek = 20.0, max_Imk = 20.0
            , extinction_rescale = false)
    for ω in ωs]

    k_eff_φs = wavenumber_low_volfrac(ωs, medium, species; tol=tol,
        hankel_order=hankel_order,
        radius_multiplier = radius_multiplier)

    wave_eff_φs = [
        EffectiveWave(ωs[i], k_eff_φs[i], medium, species; tol = tol,
            hankel_order=hankel_order, θin = θin,
            extinction_rescale = true,
            radius_multiplier = radius_multiplier)
    for i in eachindex(ωs)]

    @test maximum(abs(wave_effs_arr[i][1].k_eff/k_eff_φs[i] - 1.0) for i in eachindex(ωs)) < 5e-8

    ds = [abs(dot(wave_eff_φs[i].amplitudes[:], wave_effs_arr[i][1].amplitudes[:])) for i in eachindex(ωs)]
    ds2 = [norm(wave_eff_φs[i].amplitudes[:])*norm(wave_effs_arr[i][1].amplitudes[:]) for i in eachindex(ωs)]
    @test maximum(ds./ds2 .- 1.0) < tol

    match_ws = [
        MatchWave(ωs[i], medium, specie; θin = θin,
            radius_multiplier = radius_multiplier,
            max_size=150,
            tol = tol, wave_effs = wave_effs_arr[i][1:1])
    for i in eachindex(ωs)]

    @test maximum(match_error.(match_ws)) < 1e-6

    avg_wave_φs = [
        AverageWave(match_ws[i].average_wave.x, wave_eff_φs[i])
    for i in eachindex(ωs)]

    Rs_near = [reflection_coefficient(ωs[i], match_ws[i].average_wave, medium, specie) for i in eachindex(ωs)]
    Rs_near_φ = [reflection_coefficient(ωs[i], avg_wave_φs[i], medium, specie) for i in eachindex(ωs)]
    @test maximum(abs.(Rs_near_φ - Rs_near)) < 5e-6

    Rs = [reflection_coefficient(ωs[i], match_ws[i], medium, specie) for i in eachindex(ωs)]
    Rs_φ = [reflection_coefficient(ωs[i], wave_eff_φs[i], medium, species) for i in eachindex(ωs)]
    @test maximum(abs.(Rs_near_φ - Rs_near)) < 5e-6

    @test maximum(
        norm(match_ws[i].average_wave.amplitudes[:] .- avg_wave_φs[i].amplitudes[:])/norm(avg_wave_φs[i].amplitudes[:])
    for i in eachindex(ωs)) < 1e-3
end
