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
        MatchWave(ω, medium, species[i]; hankel_order=hankel_order,
        θin = θin, tol = tol,
        wave_effs = wave_effs_arr[i])
    for i in eachindex(species)]


    for i in eachindex(species)
        x = linspace(match_ws[i].x_match[end],200*match_ws[i].x_match[end],200)
        avg_low = AverageWave(x, wave_eff_lows[i])
        avg = AverageWave(x, match_ws[i].effective_waves)
        @test maximum(abs.(avg.amplitudes[:]./avg_low.amplitudes[:] .- 1)) < 210*tol
    end
end
