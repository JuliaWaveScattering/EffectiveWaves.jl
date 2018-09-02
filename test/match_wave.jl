@testset "match wave low volfrac" begin

## Low volume fraction
    medium = Medium(1.0,1.0+0.0im)
    specie = Specie(ρ=0.5, r=0.4, c=0.5, volfrac=0.001)
    species = [specie]
    θin = 0.2
    tol = 1e-8
    hankel_order = 2

    ωs = [0.2,1.2]
    radius_multiplier = 1.005

    wave_effs_arr = [
        effective_waves(ω, medium, [specie];
            hankel_order=hankel_order, tol = tol,  θin = θin,
            radius_multiplier = radius_multiplier, mesh_points = 10, mesh_size = 2.0 #, max_Rek = 20.0, max_Imk = 20.0
            , extinction_rescale = false)
    for ω in ωs]

    k_eff_φs = wavenumber_low_volfrac(ωs, medium, species; tol=tol,
        hankel_order=hankel_order,
        radius_multiplier = radius_multiplier)
    wave_eff_φs = [
        EffectiveWave(ωs[i], k_eff_φs[i], medium, species; tol = tol,
            hankel_order=hankel_order, θin = θin, radius_multiplier = radius_multiplier)
    for i in eachindex(ωs)]

    @test maximum(abs(wave_effs_arr[i][1].k_eff/k_eff_φs[i] - 1.0) for i in eachindex(ωs)) < 5*tol

    ds = [abs(dot(wave_eff_φs[i].amplitudes[:], wave_effs_arr[i][1].amplitudes[:])) for i in eachindex(ωs)]
    ds2 = [norm(wave_eff_φs[i].amplitudes[:])*norm(wave_effs_arr[i][1].amplitudes[:]) for i in eachindex(ωs)]
    @test maximum(ds./ds2 .- 1.0) < tol

    match_ws = [
        MatchWave(ωs[i], medium, specie; θin = θin,
            radius_multiplier = radius_multiplier,
            tol = tol, wave_effs = wave_effs_arr[i])
    for i in eachindex(ωs)]

    avg_wave_φs = [
        AverageWave(match_ws[i].average_wave.x, wave_eff_φs[i])
    for i in eachindex(ωs)]

    @test maximum(
        norm(match_ws[i].average_wave.amplitudes[:] .- avg_wave_φs[i].amplitudes[:])/norm(avg_wave_φs[i].amplitudes[:])
    for i in eachindex(ωs)) < 10*sqrt(tol)
end

@testset "match wave low frequency" begin

## Low frequency
    medium = Medium(1.0,1.0+0.0im)
    volfracs = 0.1:0.2:0.5
    species = [
        Specie(ρ=5.,r=0.002, c=v, volfrac=v)
    for v in volfracs]

    θin = -0.2
    tol = 1e-7
    hankel_order = 1

    ω = 0.002

    k_eff_lows = map(species) do s
        med_eff = effective_medium(medium, [s])
        ω/med_eff.c
    end
    wave_eff_lows = [
        EffectiveWave(ω, k_eff_lows[i], medium, [species[i]]; tol = tol,
            hankel_order=hankel_order, θin = θin)
    for i in eachindex(species)]

    wave_effs_arr = [
        effective_waves(ω, medium, [s];
            hankel_order=hankel_order, tol = tol,  θin = θin,
            extinction_rescale = false)
    for s in species]

    match_ws = [
        MatchWave(ω, medium, species[i]; hankel_order=hankel_order, θin = θin, tol = tol, wave_effs = wave_effs_arr[i])
    for i in eachindex(species)]
    @test true
    x = linspace(0.,2*match_ws[1].x_match[end],200)
    i = 3
    plot(match_ws[i])
    plot!(x, [wave_eff_lows[i]], linestyle=:dash,linewidth=2)
end
