using EffectiveWaves
using LinearAlgebra, Statistics, Test

@testset "Pair-correlation tests" begin

    particle_medium = Acoustic(3; ρ=0.0, c=0.0);

    exclusion_distance = 1.01

    s1 = Specie(
        particle_medium, Sphere(rand()+0.2), exclusion_distance = exclusion_distance
    );

    # minimal distance between particle centres
    a12 = 2.0 * s1.exclusion_distance * outer_radius(s1)

    r1s = (20 * a12) .* rand(1000);
    r2s = (20 * a12) .* rand(1000);
    cosθs = rand(-1.0:0.1:1.0,1000);

    pair_corr_inf(z) = hole_correction_pair_correlation([0.0,0.0,0.0],s1, [0.0,0.0,z],s1)
    pair_radial = EffectiveWaves.pair_radial_fun(pair_corr_inf, a12; polynomial_order = 60)

    p_rads = pair_radial.(r1s,r2s,cosθs)
    p_infs = pair_corr_inf.(sqrt.(r1s.^2 .+ r2s.^2 .- 2r1s .* r2s .* cosθs))

    maximum(abs.(p_infs))
    minimum(abs.(p_infs))

    @test mean(abs.(p_rads - p_infs)) < 1e-3
    @test maximum(abs.(p_rads - p_infs)) < 0.55 # infinite sharp function may have point wise error of 0.5, but exceptionally rare.
end
