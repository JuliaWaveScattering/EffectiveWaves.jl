@testset "match wave low frequency" begin

## Low frequency
    medium = Medium(1.0,1.0+0.0im)
    volfracs = [0.1,0.5]
    species = [
        Specie(ρ=5.,r=0.002, c=v, volfrac=v)
    for v in volfracs]

    ω = 0.002
    θin = -0.2
    tol = 1e-8
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
            hankel_order=hankel_order,
            num_wavenumbers = 1,
            tol = tol,  θin = θin,
            extinction_rescale = false)
    for s in species]

    @test norm([ ws[1].k_eff for ws in wave_effs_arr] - k_eff_lows) < tol

    match_ws = [
        MatchWave(ω, medium, species[i]; hankel_order=hankel_order,
        max_size = 60,
        θin = θin, tol = tol,
        wave_effs = wave_effs_arr[i])
    for i in eachindex(species)]

   @test maximum(match_error.(match_ws)) < tol

    map(eachindex(species)) do i
        x = linspace(match_ws[i].x_match[end],2pi/real(match_ws[i].effective_waves[1].k_eff),200)
        avg_low = AverageWave(x, wave_eff_lows[i])
        avg = AverageWave(x, match_ws[i].effective_waves)
        @test norm(avg.amplitudes[:] - avg_low.amplitudes[:]) < tol
        @test maximum(abs.(avg.amplitudes[:] - avg_low.amplitudes[:])) < tol
    end
end
