@testset "match purely numerical solution" begin

## high attenuating material
    medium = Medium(1.0,1.0+0.0im)
    # cs = [0.1,0.5]
    cs = [0.5]
    species = [
        Specie(ρ = c, c = 0.6-c, r = 2*c, volfrac=c)
    for c in cs]

    ω = 1.1
    k = ω/medium.c
    θin = 0.3
    tol = 1e-8
    hankel_order = 2
    radius_multiplier = 1.005

    wave_effs_arr = [
        effective_waves(ω, medium, [s];
            radius_multiplier=radius_multiplier,
            hankel_order=hankel_order,
            mesh_points=6,
            num_wavenumbers=28,
            tol = tol,  θin = θin,
            extinction_rescale = false)
    for s in species];

   # use only 20 least attenuating
   wave_effs_arr = [w[1:20] for w in wave_effs_arr]

    match_ws = [
        MatchWave(ω, medium, species[i];
            radius_multiplier=radius_multiplier,
            hankel_order=hankel_order,
            θin = θin, tol = tol,
            wave_effs = wave_effs_arr[i],
            max_size=80)
            # x = 0.0:(radius_multiplier*2*species[i].r/100):4)
    for i in eachindex(species)];

    @test maximum(match_error.(match_ws)) < 1e-4

    avgs = [
        AverageWave(ω, medium, species[i];
                radius_multiplier=radius_multiplier,
                hankel_order=hankel_order,
                tol = tol, θin = θin,
                wave_effs = wave_effs_arr[i], max_size=700)
    for i in eachindex(species)]

    R_ms = [reflection_coefficient(ω, match_ws[i], medium, species[i]) for i in eachindex(species)]
    R_ds = [reflection_coefficient(ω, avgs[i], medium, species[i]) for i in eachindex(species)]

    avg_eff = AverageWave(match_ws[1].x_match[end]:0.002:40, match_ws[1].effective_waves);
    R1 = reflection_coefficient(ω, avg_eff, medium, species[1])
    R2 = reflection_coefficient(ω, match_ws[1].effective_waves, medium, species; x=avg_eff.x[1])
    @test norm(R1 - R2) < 1e-7

    @test norm(R_ms - R_ds) < 5e-4

    map(eachindex(species)) do i
        j0 = findmin(abs.(avgs[i].x .- match_ws[i].x_match[1]))[2]
        x0 = avgs[i].x[j0+1:end]
        avg_m = AverageWave(x0, match_ws[i].effective_waves)
        maximum(abs.(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:]))
        @test norm(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:])/norm(avg_m.amplitudes[:]) < 4e-3
        @test maximum(abs.(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:])) < 2e-3
    end
end
