using EffectiveWaves
using LinearAlgebra, Statistics, Test

@testset "Pair-correlation tests" begin

    particle_medium = Acoustic(3; ρ=0.0, c=0.0);

    exclusion_distance = 1.01

    s1 = Specie(
        particle_medium, Sphere(0.8), exclusion_distance = exclusion_distance
    );

    # minimal distance between particle centres
    a12 = 2.0 * s1.exclusion_distance * outer_radius(s1)
    R = 10.0
    polynomial_order = 40

    pair_corr_inf(z) = hole_correction_pair_correlation([0.0,0.0,0.0],s1, [0.0,0.0,z],s1)

    # import EffectiveWaves: smooth_pair_corr_distance

    pair_corr_inf_smooth = smooth_pair_corr_distance(pair_corr_inf, a12; smoothing = 1.0, max_distance = 2R, polynomial_order = polynomial_order)

    # using Plots
    # rs = 0.0:0.1:(2R)
    # plot(rs,pair_corr_inf_smooth.(rs))
    # plot!(rs,pair_corr_inf.(rs))

    pair_radial = pair_radial_fun(pair_corr_inf_smooth, a12;
        sigma_approximation = false,
        polynomial_order = polynomial_order
    )

    r1 = 4.0; r2 = 4.2;
    cosθs = LinRange(0.6,1.0,200);
    p_rads = pair_radial.(r1,r2,cosθs);
    p_infs = pair_corr_inf_smooth.(sqrt.(r1^2 .+ r2.^2 .- 2r1 .* r2 .* cosθs))

    maximum(abs.(p_rads - p_infs))

    # plot(cosθs, p_rads)
    # plot!(cosθs,p_infs)

    @test maximum(abs.(p_rads - p_infs)) < 1e-3

    polynomial_order = 90;
    pair_corr_inf(z) = hole_correction_pair_correlation([0.0,0.0,0.0],s1, [0.0,0.0,z],s1) #* (1 + sin(5*abs(z-a12)) * exp(-z))

    pair_radial = pair_radial_fun(pair_corr_inf, a12;
        sigma_approximation = true,
        # sigma = true,
        polynomial_order = polynomial_order,
        mesh_size = 12polynomial_order + 1
    )

    # polynomial_order = 10
    # sigmas = [one(T); sin.(pi .* ls[2:end] ./ (polynomial_order+1)) ./ (pi .* ls[2:end] ./ (polynomial_order+1))]
    # plot(sigmas)

    p_rads = pair_radial.(r1,r2,cosθs);
    p_infs = pair_corr_inf.(sqrt.(r1^2 .+ r2.^2 .- 2r1 .* r2 .* cosθs))
    # using Plots
    #
    # plot(cosθs,p_rads)
    # plot!(cosθs,p_infs)

    r1s = (20 * a12) .* rand(1000);
    r2s = (20 * a12) .* rand(1000);
    cosθs = LinRange(-1.0,1.0,1000);

    p_rads = pair_radial.(r1s,r2s,cosθs);
    p_infs = pair_corr_inf.(sqrt.(r1s.^2 .+ r2s.^2 .- 2r1s .* r2s .* cosθs));

    @test mean(abs.(p_rads - p_infs)) < 1e-3
    @test maximum(abs.(p_rads - p_infs)) < 0.52 # infinite sharp function may have point wise error of 0.5, but exceptionally rare.
end
