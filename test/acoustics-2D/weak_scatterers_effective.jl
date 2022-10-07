@testset "weak scatterers" begin

    # angular frequencies
    ωs = [0.001,20.]
    # ωs = [20.] # use only one frequency to speed up test

    # tolerance
    tol = 1e-5;
    basis_order = 2;

    # background medium
    medium = Acoustic(2; ρ=1.0, c=1.0)
    ms = MultipleScattering

    p1 = Particle(Acoustic(2; ρ=10.0, c=12.),ms.Circle(0.01))
    p2 = Particle(Acoustic(2; ρ=3.0, c=2.0),ms.Circle(0.2))

    species = [
        Specie(p1; volume_fraction=0.03),
        Specie(p2; volume_fraction=0.04)
    ]

    # wavenumbers
    eff_medium = effective_medium(medium, species)

    k_eff_lows = ωs./eff_medium.c
    k_eff_φs = wavenumber_low_volumefraction(ωs, medium, species; basis_order = basis_order)

    k_effs_arr = [
        wavenumbers(ω, medium, species;
            num_wavenumbers=1, tol = tol, basis_order = basis_order)
    for ω in ωs]

    inds = [argmin(abs.(k_effs_arr[i] .- k_eff_φs[i])) for i in eachindex(ωs)]
    k_effs = [k_effs_arr[i][inds[i]] for i in eachindex(inds)]

    # k_effs =
    # 0.0010101407549128578 + 7.824113157942236e-13im
    # 20.207827596241156 + 0.11344062283733775im

    @test norm( (k_effs - k_eff_φs) ./ norm.(k_eff_φs) ) < 0.001
    @test norm( (k_effs - k_eff_lows) ./ norm.(k_eff_lows) ) < 0.01
    @test norm(k_effs[1] - k_eff_lows[1])/norm(k_effs[1]) < 1e-7

    normal = [-1.0,0.0] # an outward normal to the surface
    material = Material(Halfspace(normal),species)

    # define a plane wave source travelling at a 45 degree angle in relation to the material
    source = PlaneSource(medium, -normal)

    wave_effs_φs = [WaveMode(ωs[i], k_eff_φs[i], source, material; tol = tol) for i in eachindex(ωs)]

    # reflection coefficient
    Rs2 = reflection_coefficients(ωs, source, material; tol = tol, num_wavenumbers = 1, basis_order = basis_order)

    # if the wavenumber k_effs are already calculated, the below is faster
    # pairs each ω in ωs with each wave in wave_effs2 to calculate the refleciton coefficients
    Rs = map(eachindex(ωs)) do i
        w_effs = [
            WaveMode(ωs[i], k_eff, source, material; tol = tol, basis_order=basis_order)
        for k_eff in k_effs_arr[i]]
        reflection_coefficient(ωs[i], w_effs, source, material; tol = tol)
    end

    @test Rs ≈ Rs2 # same keywords should lead to the same reflection coefficient, except for numerically rounding

    Rs1 = map(eachindex(ωs)) do i
        w_eff = WaveMode(ωs[i], k_effs[i][1], source, material; tol = 1e-9)
        reflection_coefficient(ωs[i], w_eff, source, material)
    end

    # Direct incidence
    R_low = reflection_transmission_coefficients(ωs[1], source, eff_medium, material.shape)[1]
    Rs_φs = reflection_coefficients(ωs, wave_effs_φs, source, material; tol=1e-9)
    # the below takes a low-volfrac expansion for both the wavenumber and reflection coefficient
    Rs_φs2 = reflection_coefficient_low_volumefraction(ωs, source, material; basis_order = basis_order)

    @test maximum(abs.(Rs_φs - Rs1)) < 1e-4 # already relative to incident wave amplitude = 1
    @test maximum(abs.(Rs_φs2 - Rs1)) < 1e-3
    @test abs(R_low - Rs1[1]) < 1e-7

    # Vary angle of incidence θin. Beware of extreme angles θ ~= -pi/2 or pi/2.
    θs = (-π/2.1):0.3:(π/2.1)
    sources = [PlaneSource(medium, [cos(θ),sin(θ)]) for θ in θs]

    R_low = map(sources) do s
        reflection_transmission_coefficients(ωs[1], s, eff_medium, material.shape)[1]
    end

    Rs = map(sources) do s
        wave_effs_2 = [
            WaveMode(ωs[i], k_effs[i], s, material; tol = tol, basis_order = basis_order)
        for i in eachindex(ωs)]
        reflection_coefficients(ωs, wave_effs_2, s, material)
    end
    Rs_φs = [
        reflection_coefficient_low_volumefraction(ωs, s, material; basis_order = basis_order)
    for s in sources]

    len = length(ωs)

    @test maximum(abs(R_low[i] - Rs[i][1]) for i in eachindex(R_low)) < 1e-7
    @test maximum(norm(R)/len for R in (Rs_φs - Rs)[2:end-2]) < 5e-3
end
