using EffectiveWaves, Test
using LinearAlgebra

# this is for high volume and frequency
@testset "compare wienger-hopf and match method average scattering coefficient" begin

    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    specie = Specie(Particle(
        Acoustic(2; ρ=0.0, c=0.0), ms.Circle(0.5));
        volume_fraction=0.3
    )

    ω = 1.0

    tol = 1e-8
    basis_order=0
    θ = pi/4

    k_effs = wavenumbers(ω, medium, [specie];
        tol=tol, basis_order = basis_order,
        num_wavenumbers = 50);

    wave_effs = [
        EffectivePlaneWaveMode(ω, k_eff, medium, [specie];
            basis_order = basis_order,
            tol = tol, extinction_rescale=false,
            method = :WienerHopf,
            θin=θ
        )
    for k_eff in k_effs];

    match_ws = MatchPlaneWaveMode(ω, medium, specie;
        basis_order = basis_order,
        tol = tol, wave_effs = wave_effs[1:10],
        θin=θ,
        max_size = 800,
    );

    inds = Int.(round.(LinRange(1,length(wave_effs),5)))
    errors = map(inds) do i
        avg_WH = DiscretePlaneWaveMode(match_ws.discrete_wave.x[1:200], wave_effs[1:i]);
        norm(avg_WH.amplitudes - match_ws.discrete_wave.amplitudes[1:200,:,:])/norm(match_ws.discrete_wave.amplitudes[1:200,:,:])
    end

    # errors should be monotonically decreasing
    @test sort(errors; rev=true) == errors
    @test errors[end] < 0.05
end
