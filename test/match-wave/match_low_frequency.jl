@testset "match wave low frequency" begin

## Low frequency
    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    volfracs = [0.1,0.5]
    species = [
        Specie(Particle(
            Acoustic(2; ρ=5., c=v), ms.Circle(0.002));
            volume_fraction=v
        )
    for v in volfracs]

    ω = 0.002
    θin = -0.2
    tol = 1e-7
    basis_order = 1

    k_eff_lows = map(species) do s
        med_eff = effective_medium(medium, [s])
        ω/med_eff.c
    end

    normal = [-1.0,0.0] # an outward normal to the surface
    materials = [Material(Halfspace(normal),s) for s in species]
    source = PlaneSource(medium, [cos(θin),sin(θin)])

    wave_eff_lows = [
        WaveMode(ω, k_eff_lows[i], source, materials[i]; tol = tol,
            basis_order=basis_order)
    for i in eachindex(species)]

    wave_effs_arr = [
        WaveModes(ω, source, materials[1];
            basis_order = basis_order,
            num_wavenumbers = 1,
            tol = tol,
            extinction_rescale = false),
        WaveModes(ω, source, materials[2];
            basis_order = basis_order,
            num_wavenumbers = 1,
            tol = tol,
            extinction_rescale = false)
    ]
    # for i in eachindex(species)]
    # Using the error causes a segmentation fault

    @test norm([ws[1].wavenumber for ws in wave_effs_arr] - k_eff_lows) < tol

    match_ws = [
        MatchPlaneWaveMode(ω, source, materials[i]; basis_order=basis_order,
            max_size = 60,
            tol = tol,
            wave_effs = wave_effs_arr[i])
    for i in eachindex(species)]

   @test maximum(match_error(match_ws[i],materials[i].shape) for i in eachindex(species)) < tol

    map(eachindex(species)) do i
        sx = abs(real(match_ws[i].PlaneWaveModes[1].wavenumber))
        x = LinRange(match_ws[i].x_match[end],2pi / sx,200)
        avg_low = DiscretePlaneWaveMode(x, wave_eff_lows[i], materials[i].shape)
        avg = DiscretePlaneWaveMode(x, match_ws[i].PlaneWaveModes, materials[i].shape)
        @test norm(avg.amplitudes[:] - avg_low.amplitudes[:]) < tol
        @test maximum(abs.(avg.amplitudes[:] - avg_low.amplitudes[:])) < tol
    end
end
