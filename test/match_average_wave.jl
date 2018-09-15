@testset "match purely numerical solution" begin

## high attenuating material
    medium = Medium(1.0,1.0+0.0im)
    # cs = [0.1,0.5]
    cs = [0.5]
    species = [
        Specie(ρ = c, c = 0.6-c, r = 2*c, volfrac=c)
    for c in cs]

    ω = 1.1
    θin = 0.3
    tol = 1e-8
    hankel_order = 2

    wave_effs_arr = [
        effective_waves(ω, medium, [s];
            hankel_order=hankel_order,
            mesh_points=6,
            num_wavenumbers=10,
            tol = tol,  θin = θin,
            extinction_rescale = false)
    for s in species];

    match_ws = [
        MatchWave(ω, medium, species[i];
            hankel_order=hankel_order,
            θin = θin, tol = tol,
            wave_effs = wave_effs_arr[i],
            max_size=80)
    for i in eachindex(species)]

    avgs = [
        AverageWave(ω, medium, species[i];
                hankel_order=hankel_order,
                tol = tol, θin = θin,
                wave_effs = wave_effs_arr[i], max_size=500)
    for i in eachindex(species)]

    map(eachindex(species)) do i
        j0 = findmin(abs.(avgs[i].x .- match_ws[i].x_match[1]))[2]
        x0 = avgs[i].x[j0+1:end]
        avg_m = AverageWave(x0, match_ws[i].effective_waves)
        maximum(abs.(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:]))
        @test norm(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:])/norm(avg_m.amplitudes[:]) < 4e-3
        @test maximum(abs.(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:])) < 4e-3
    end
end
