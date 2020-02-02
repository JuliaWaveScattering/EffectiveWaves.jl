using EffectiveWaves, Test
using LinearAlgebra

# this is for high volume and frequency
@testset "compare wienger-hopf and match method average scattering coefficient" begin

    medium = Acoustic(2; ρ=1.0, c=1.0)
    ω = 1.0

    tol = 1e-9
    basis_order=0
    θ = pi/4

    ms = MultipleScattering
    specie = Specie(Particle(
        Acoustic(2; ρ=0.0, c=0.0), ms.Circle(0.5));
        volume_fraction=0.3
    )

    normal = [-1.0,0.0]; # an outward normal to the surface
    material = Material(Halfspace(normal),specie);
    source = PlaneSource(medium, [cos(θ),sin(θ)]);

    k_effs = wavenumbers(ω, medium, [specie];
        tol=tol, basis_order = basis_order,
        num_wavenumbers = 20);

    wave_effs = [
        effective_wavemode(ω, k_eff, source, material;
            basis_order = basis_order,
            tol = tol, extinction_rescale=false,
            method = :WienerHopf
        )
    for k_eff in k_effs];

    match_ws = MatchPlaneWaveMode(ω, source, material;
        basis_order = basis_order,
        tol = tol/10., wave_effs = wave_effs[1:15],
        max_size = 1000,
    );

    inds = Int.(round.(LinRange(1,length(wave_effs),5)))
    errors = map(inds) do i
        avg_WH = DiscretePlaneWaveMode(match_ws.discrete_wave.x[1:200], wave_effs[1:i], material.shape);
        norm(avg_WH.amplitudes - match_ws.discrete_wave.amplitudes[1:200,:,:])/norm(match_ws.discrete_wave.amplitudes[1:200,:,:])
    end

    # errors should be monotonically decreasing
    @test sort(errors; rev=true) == errors
    @test errors[end] < 0.07 # previously was < 0.05 before when transmission angle was limited to propagate into the material
end
