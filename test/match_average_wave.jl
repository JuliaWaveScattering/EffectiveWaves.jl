using EffectiveWaves, Test

@testset "match purely numerical solution" begin

## high attenuating material
    medium = Acoustic(2; ρ=1.0, c=1.0)
    # cs = [0.1,0.5]
    cs = [0.1,0.5]
    ms = MultipleScattering

    species = [
        Specie(Particle(
            Acoustic(2; ρ=c, c=0.6-c), ms.Circle(2.0*c));
            volume_fraction=c
        )
    for c in cs]

    ω = 1.1
    k = ω/medium.c
    θin = 0.3
    tol = 1e-7
    basis_order = 2

    wave_effs_arr = [
        effective_waves(ω, medium, [s];
            basis_order=basis_order,
            mesh_points=5,
            num_wavenumbers=5,
            tol = tol,  θin = θin,
            extinction_rescale = false)
    for s in species];

   # use only 5 least attenuating
   wave_effs_arr = [w[1:5] for w in wave_effs_arr]

    match_ws = [
        MatchWave(ω, medium, species[i];
            basis_order=basis_order,
            θin = θin, tol = tol,
            wave_effs = wave_effs_arr[i],
            max_size=80)
    for i in eachindex(species)];

    @test maximum(match_error.(match_ws)) < 1e-4

    avgs = [
        AverageWave(ω, medium, species[i];
                basis_order=basis_order,
                tol = tol, θin = θin,
                wave_effs = wave_effs_arr[i], max_size=700)
    for i in eachindex(species)]

    R_ms = [reflection_coefficient(ω, match_ws[i], medium, species[i]) for i in eachindex(species)]
    R_ds = [reflection_coefficient(ω, avgs[i], medium, species[i]) for i in eachindex(species)]

    avg_eff = AverageWave(match_ws[2].x_match[end]:0.002:40, match_ws[2].effective_waves);
    R1 = reflection_coefficient(ω, avg_eff, medium, species[2])
    R2 = reflection_coefficient(ω, match_ws[2].effective_waves, medium, [species[2]]; x=avg_eff.x[1])
    @test norm(R1 - R2) < 1e-7

    @test maximum(abs.(R_ms - R_ds)) < 8e-4

    map(eachindex(species)) do i
        j0 = findmin(abs.(avgs[i].x .- match_ws[i].x_match[1]))[2]
        x0 = avgs[i].x[j0+1:end]
        avg_m = AverageWave(x0, match_ws[i].effective_waves)
        maximum(abs.(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:]))
        @test norm(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:])/norm(avg_m.amplitudes[:]) < 6e-3
        @test maximum(abs.(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:])) < 2e-3
    end
end
