using EffectiveWaves, Test
using LinearAlgebra

# this is for high volume and frequency
@testset "compare wienger-hopf and match method average scattering coefficient" begin

    medium = Medium(1.0,1.0+0.0im)

    ω = 1.0

    specie = Specie(ρ=0.0, r=0.5, c=0.0, volfrac = 0.30)
    radius_multiplier = 1.001

    tol = 1e-8
    hankel_order=0
    θ = pi/4

    k_effs = wavenumbers(ω, medium, [specie];
        radius_multiplier = radius_multiplier,
        tol=tol, hankel_order = hankel_order,
        num_wavenumbers = 50);

    wave_effs = [
        EffectiveWave(ω, k_eff, medium, [specie];
            hankel_order = hankel_order,
            radius_multiplier = radius_multiplier,
            tol = tol, extinction_rescale=false,
            method = :WienerHopf,
            θin=θ
        )
    for k_eff in k_effs];

    match_ws = MatchWave(ω, medium, specie;
        radius_multiplier = radius_multiplier,
        hankel_order = hankel_order,
        tol = tol, wave_effs = wave_effs[1:10],
        θin=θ,
        max_size = 800,
    );

    inds = Int.(round.(LinRange(1,length(wave_effs),5)))
    errors = map(inds) do i
        avg_WH = AverageWave(match_ws.average_wave.x[1:200], wave_effs[1:i]);
        norm(avg_WH.amplitudes - match_ws.average_wave.amplitudes[1:200,:,:])/norm(match_ws.average_wave.amplitudes[1:200,:,:])
    end

    # errors should be monotonically decreasing
    @test sort(errors; rev=true) == errors
    @test errors[end] < 0.05
end
