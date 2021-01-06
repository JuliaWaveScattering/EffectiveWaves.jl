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
    tol = 1e-8
    basis_order = 2

    normal = [-1.0,0.0] # an outward normal to the surface
    materials = [Material(Halfspace(normal),s) for s in species]
    source = PlaneSource(medium, [cos(θin),sin(θin)])

    material = materials[1]

    # import StaticArrays: SVector
    #
    # function unreachable_test(ω::T, source::AbstractSource{T}, material::Material{Dim,S,Sps}; kws...) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}}
    #
    #     k_effs = rand(2) + rand(2) .* im
    #     wave_effs = map(k_effs) do k_eff
    #         # amp = Complex{Float64}[0.000326919-0.000330135im; -0.000605635-0.0396968im; -0.880755+0.456138im; -0.000605635-0.0396968im; 0.000326919-0.000330135im];
    #         amp = rand(5) + rand(5) .* im;
    #         # wavevector = Complex{Float64}[1.62185+1.33695im, 0.325072+0.0im]
    #         wavevector = rand(2) + rand(2) .* im
    #         basis_order = 2
    #         EffectivePlaneWaveMode(ω, basis_order, SVector(wavevector...), amp)
    #         # WaveMode(ω, k_eff, source, material; kws...)
    #     end
    #
    #     return wave_effs
    # end
    #
    # data = [
    #     unreachable_test(ω, source, materials[i])
    # for i in eachindex(species)]

    wave_effs_arr = Vector{Vector{EffectivePlaneWaveMode{Float64,2}}}(undef,length(species))

    i = 1
    wave_effs_arr[i] = WaveModes(ω, source, materials[i];
        basis_order=basis_order,
        box_k = [[-18.0,18.0],[0.0,20.0]],
        num_wavenumbers = 6,
        tol = tol)

    # k_effs = wavenumbers_bisection(ω, medium, [species[1]];
    #     basis_order = basis_order,
    #     box_k = [[-18.0,18.0],[0.0,20.0]],
    #     num_wavenumbers = 8,
    #     tol = tol)

    i = 2
    wave_effs_arr[i] = WaveModes(ω, source, materials[i];
        basis_order=basis_order,
        num_wavenumbers=6,
        tol = tol);
        # extinction_rescale = false)

    # Causes unreachable error..
    # for i in eachindex(species)
    #     wave_effs_arr[i] = WaveModes(ω, source, materials[i];
    #         basis_order=basis_order,
    #         mesh_points=5,
    #         num_wavenumbers=5,
    #         tol = tol,
    #         extinction_rescale = false)
    # end

   # use only 6 least attenuating
   wave_effs_arr2 = [w[1:6] for w in wave_effs_arr];

    match_ws = [
        MatchPlaneWaveMode(ω, source, materials[i];
            basis_order=basis_order,
            tol = tol,
            wave_effs = wave_effs_arr2[i],
            max_size=300)
            # max_size=480)
    for i in eachindex(species)];

    @test maximum(match_error(match_ws[i],materials[i].shape) for i in eachindex(species)) < 200*tol

    avgs = [
        DiscretePlaneWaveMode(ω, source, materials[i];
                basis_order=basis_order,
                tol = tol, max_size=600,
                wave_effs = wave_effs_arr2[i])
    for i in eachindex(species)]

    R_ms = [reflection_coefficient(ω, match_ws[i], source, materials[i]) for i in eachindex(species)]
    R_ds = [reflection_coefficient(ω, avgs[i], source, materials[i]) for i in eachindex(species)]
    @test maximum(abs.(R_ms - R_ds)) < 5e-4

    avg_eff = DiscretePlaneWaveMode(match_ws[2].x_match[end]:0.002:40, match_ws[2].PlaneWaveModes, material.shape);
    R1 = reflection_coefficient(ω, avg_eff, source, materials[2])
    R2 = reflection_coefficient(ω, match_ws[2].PlaneWaveModes, source, materials[2]; x=avg_eff.x[1])
    @test norm(R1 - R2) < 1e-7

    map(eachindex(species)) do i
        j0 = findmin(abs.(avgs[i].x .- match_ws[i].x_match[1]))[2]
        x0 = avgs[i].x[j0+1:end]
        avg_m = DiscretePlaneWaveMode(x0, match_ws[i].PlaneWaveModes, materials[i].shape)
        maximum(abs.(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:]))
        @test norm(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:])/norm(avg_m.amplitudes[:]) < 1e-2
        @test maximum(abs.(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:])) < 3e-4
        maximum(abs.(avgs[i].amplitudes[j0+1:end,:,:][:] - avg_m.amplitudes[:]))
    end
end
