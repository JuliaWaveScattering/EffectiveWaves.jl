using EffectiveWaves, ClassicalOrthogonalPolynomials
using LinearAlgebra, Statistics, Test

@testset "Discrete radial system for sphere scattering" begin

## Set parameters
    particle_medium = Acoustic(3; ρ=0.1, c=0.1);
    medium = Acoustic(3; ρ=1.0, c=1.0);

    R = 5.0; r = 1.0

    separation_ratio = 1.02
    kas = [0.01,0.2]
    ks = kas ./ r

    vol_fraction = 0.15

    basis_orders = Int.(round.(4.0 .* kas)) .+ 1
    basis_field_orders = Int.(round.(3.0 .* ks .* R)) .+ 1

    ωs = ks .* real(medium.c)

    s1 = Specie(
        particle_medium, Sphere(r);
        number_density = vol_fraction / volume(Sphere(r)),
        exclusion_distance = separation_ratio
    );

    species = [s1]
    # species = [s1,s1]

    region_shape = Sphere([0.0,0.0,0.0], R)
    material = Material(Sphere(R),species);

    sourceradial = regular_spherical_source(medium, [1.0+0.0im];
       position = [0.0,0.0,0.0], symmetry = RadialSymmetry{3}()
    );

    a12 = 2.0 * s1.exclusion_distance * outer_radius(s1)
    rs = 0.0:0.02:(R - outer_radius(s1));
    xs = [
        radial_to_cartesian_coordinates([r,0.5,0.6])
    for r in rs];



## The full discrete system does not calculate the translation matrix for |r1 - r2| <= a12 , whereas the radial method does. To get an exact match between the methods we will use a special pair-correlation which does not allow two particles on the same spherical annulus

    polynomial_order = 1
    P = Legendre()

    function gls_radial_simple(r1,r2)
        if abs(r1-r2) < a12
            zeros(polynomial_order+1)
        else
            1 ./ (1:(polynomial_order+1)) .* cos(abs(r1-r2) * pi/(2R))
        end
    end

    function pair_radial_simple(r1,r2,u)
        Pus = P[u, 1:(polynomial_order + 1)] .* (2 .* (0:polynomial_order) .+ 1) ./ (4pi)
        return sum(Pus .* gls_radial_simple(r1,r2))
    end

    function pair_corr_simple(x1,s1,x2,s2)
        if norm(x1) < 1e-12 || norm(x2) < 1e-12
            pair_radial_simple(norm(x1),norm(x2),0.0)
        else
            pair_radial_simple(norm(x1),norm(x2),dot(x1,x2) / (norm(x1)*norm(x2)))
        end
    end

    discrete_field_radials = [
        EffectiveWaves.discrete_system_radial(ωs[i], sourceradial, material, Symmetry(sourceradial,material);
            basis_order = basis_orders[i],
            basis_field_order = basis_field_orders[i],
            polynomial_order = polynomial_order,
            gls_pair_radial = gls_radial_simple,
            sigma_approximation = false
        )
    for i in eachindex(ωs)];

    discrete_rad_scats = [d.coefficient_field.(xs) for d in discrete_field_radials];

    mat_dcoefs_radial = [
        material_scattering_coefficients(discrete_field_radials[i]; rtol = 1e-3,maxevals = Int(5e3))
    for i in eachindex(ωs)];

    # fully numerical method
    tmp = discrete_system(ωs[1], sourceradial, material;
        basis_order = 0, basis_field_order = 0, rtol = 1.0, maxevals = 4
    );

    # import EffectiveWaves: discrete_system
    ST = typeof(tmp);
    discrete_fields = ST[
        discrete_system(ωs[i], sourceradial, material;
            basis_order = basis_orders[i],
            basis_field_order = basis_field_orders[i],
            rtol = 5e-3, maxevals = Int(1e5),
            pair_corr = pair_corr_simple
            # pair_corr = pair_corr_smooth
        )
    for i in eachindex(ωs)];

    discrete_scats = [d.coefficient_field.(xs) for d in discrete_fields];

    mat_dcoefs = [
        material_scattering_coefficients(discrete_fields[i]; rtol = 1e-3,maxevals = Int(5e3))
    for i in eachindex(ωs)];


    errors = [
        abs(mat_dcoefs_radial[i][1] - mat_dcoefs[i][1]) / abs(mat_dcoefs[i][1])
    for i in eachindex(mat_dcoefs_radial)]

    @test errors[1] < 1e-4
    @test errors[2] < 2e-3

    errors = [
        norm(discrete_rad_scats[i] - discrete_scats[i]) / norm(discrete_rad_scats[i])
    for i in eachindex(discrete_rad_scats)];

    @test errors[1] < 1e-3
    @test errors[2] < 3e-3


## For a pair correlation that depends only on inter particle distances, the two discrete methods are only approximately the same when using a high order approximation for the pair correlation.

    # use 0 basis order to speed up code. Large order was tested above
    basis_orders = basis_orders .* 0

    polynomial_order = 60
    polynomial_order = 20
    pair_corr_inf(z) = hole_correction_pair_correlation([0.0,0.0,0.0],s1, [0.0,0.0,z],s1)

    pair_corr_inf_smooth = smooth_pair_corr_distance(
        pair_corr_inf, a12;
        smoothing = 0.0, max_distance = 2R,
        polynomial_order = polynomial_order
    )

    # using Plots
    # zs = 0.0:0.01:(2R)
    # plot(pair_corr_inf_smooth,zs)
    # plot!(pair_corr_inf,zs, linestyle=:dash)
    #
    gls_radial = gls_pair_radial_fun(pair_corr_inf, a12;
        sigma_approximation = false,
        polynomial_order = polynomial_order
    )

    # plot(abs.(gls_radial(3.0,3.0)), ylims = (0.0,2.5))

    P = Legendre{Float64}()
    function pair_radial_smooth2(r1,r2,u)
        Pus = P[u, 1:(polynomial_order + 1)] .* (2 .* (0:polynomial_order) .+ 1) ./ (4pi)

        return sum(Pus .* gls_radial(r1,r2))
    end
    #
    # plot(zs,pair_radial_smooth2.(zs,0.0,0.0))
    # plot!(zs,pair_corr_inf_smooth.(zs), linestyle=:dash)

    function pair_corr_smooth(x1,s1,x2,s2)
        if norm(x1) < 1e-12 || norm(x2) < 1e-12
            pair_radial_smooth2(norm(x1),norm(x2),0.0)
        else
            pair_radial_smooth2(norm(x1),norm(x2),dot(x1,x2) / (norm(x1)*norm(x2)))
        end
    end


## Reduced radial method
    # import EffectiveWaves: discrete_system_radial
    discrete_field_radials = [
        discrete_system_radial(ωs[i], sourceradial, material, Symmetry(sourceradial,material);
            # h12 = a12,
            basis_order = basis_orders[i],
            basis_field_order = basis_field_orders[i],
            polynomial_order = polynomial_order,
            gls_pair_radial = gls_radial,
            sigma_approximation = false
        )
    for i in eachindex(ωs)]

    discrete_rad_scats = [
        d.coefficient_field.(xs)
    for d in discrete_field_radials];

    mat_dcoefs_radial = [
        material_scattering_coefficients(discrete_field_radials[i];
            rtol = 1e-3, maxevals = Int(5e3))
    for i in eachindex(ωs)];

## fully numerical method
    tmp = discrete_system(ωs[1], sourceradial, material;
        basis_order = 0, basis_field_order = 0, rtol = 1.0, maxevals = 4
    );

    # import EffectiveWaves: discrete_system
    ST = typeof(tmp);
    discrete_fields = ST[
        discrete_system(ωs[i], sourceradial, material;
            basis_order = basis_orders[i],
            basis_field_order = basis_field_orders[i],
            rtol = 5e-3, maxevals = Int(1e5),
            pair_corr = pair_corr_smooth
            # pair_corr = hole_correction_pair_correlation
        )
    for i in eachindex(ωs)];

    discrete_scats = [
        d.coefficient_field.(xs)
    for d in discrete_fields];

    mat_dcoefs = [
        material_scattering_coefficients(discrete_fields[i]; rtol = 1e-3,maxevals = Int(5e3))
    for i in eachindex(ωs)];

    errors = [
        abs(mat_dcoefs_radial[i][1] - mat_dcoefs[i][1]) / abs(mat_dcoefs[i][1])
    for i in eachindex(mat_dcoefs_radial)]

    @test errors[1] < 1e-3
    @test errors[2] < 1e-2

    errors = [
        norm(discrete_rad_scats[i] - discrete_scats[i]) / norm(discrete_rad_scats[i])
    for i in eachindex(discrete_rad_scats)];

    # i = 2
    # fun = abs
    # plot(rs, abs.([d[1] for d in discrete_scats[i]]))
    # plot!(rs, abs.([d[1] for d in discrete_rad_scats[i]]), linestyle=:dash)

    @test errors[1] < 2e-3
    @test errors[2] < 2e-2


   #  # mat_coefs_radial = material_scattering_coefficients.(wavemodes_radial);
   # abs(mat_dcoefs_radial[i][1] - mat_dcoefs[i][1]) / abs(mat_dcoefs[i][1])
   # # abs(mat_dcoefs[i][1] - mat_coefs_radial[i][1]) / abs(mat_coefs_radial[i][1])
   # # abs(mat_dcoefs_radial[i][1] - mat_coefs_radial[i][1]) / abs(mat_coefs_radial[i][1])
   #
   # j = 1;
   # i = 2;
   # # eff_rad0s = [s[j] for s in eff_scats_radial[i]];
   #
   # df_rad0s = [s[j] for s in discrete_rad_scats[i]];
   # # df_rad1s = [s[3] for s in discrete_rad_scats[2]];
   # df0s = [d[j] for d in discrete_scats[i]];
   # # df1s = [d[3] for d in discrete_scats[2]];
   #
   # # rad_scats2 = scattered_field.(xs);
   # # df_rad2s = [s[j] for s in rad_scats2];
   #
   # zs = [x[3] for x in xs];
   # zs - norm.(xs);
   #
   # using Plots
   # gr(linewidth=2.0)
   # fun = abs;
   # plot(zs,fun.(df0s), lab = "$fun fully discrete")
   # plot!(zs,fun.(df_rad0s), linestyle = :dash, lab = "$fun discrete radial")
   # plot!(zs,fun.(eff_rad0s),
   #      linestyle = :dot,
   #      lab = "$fun eff. wave method",
   #      xlab = "radial distance",
   #      ylab = ""
   # )
end


@testset "Effective sphere scattering" begin

## Set parameters
    particle_medium = Acoustic(3; ρ=10.0, c=10.0);
    medium = Acoustic(3; ρ=1.0, c=1.0);

    R = 5.0
    r = 1.0

    separation_ratio = 1.02

    kas = [0.01,0.2]
    ks = kas ./ r

    vol_fraction = 0.12

    basis_orders = Int.(round.(4. .* kas)) .+ 1
    basis_orders = min.(basis_orders,1)

    basis_field_orders = Int.(round.(4.0 .* ks .* R)) .+ 1
    basis_field_orders = max.(basis_field_orders,2)

    ωs = ks .* real(medium.c)

    s1 = Specie(
        particle_medium, Sphere(r);
        number_density = vol_fraction / volume(Sphere(r)),
        exclusion_distance = separation_ratio
    );

    species = [s1]
    # species = [s1,s1]

    a12 = 2.0 * outer_radius(s1) * s1.exclusion_distance


## define sources and material

    psource = PlaneSource(medium, [0.0,0.0,1.0]);
    source = plane_source(medium; direction = [0.0,0.0,1.0])

    region_shape = Sphere([0.0,0.0,0.0], R)
    material = Material(Sphere(R),species);

    eff_medium = effective_medium(medium, species; numberofparticles = material.numberofparticles)
    ks_low = ωs ./ eff_medium.c

    keff_arr = [
        wavenumbers(ωs[i], medium, species;
            # num_wavenumbers = 4,
            basis_order = basis_orders[i],
            tol = 1e-7,
            numberofparticles = material.numberofparticles
        )
    for i in eachindex(ωs)]

    keffs = [ks[1] for ks in keff_arr]

## Plane wave reflection from a sphere

    ## effective waves solution
    pwavemodes_azi = [
        WaveMode(ωs[i], keffs[i], psource, material;
           basis_order = basis_orders[i],
           basis_field_order = basis_field_orders[i]
           , source_basis_field_order = basis_field_orders[i]
        )
    for i in eachindex(ωs)];

    pscat_fields = scattering_field.(pwavemodes_azi);

    ## discrete numerical solution of the average integral equations

    # increasing these parameters does lead to more accurate solutions, but convergences is slow. To increase accuracy you need to increase basis_field_order and maxevals.
    rtol = 1e-2; maxevals = Int(5e3);

    # this below is just to get the typeof tmp, to then avoid a weird unionall error which should hopefully go away when updating Julia at some point
    tmp = discrete_system(ωs[1], psource, material;
        basis_order = 0, basis_field_order = 0, rtol = 1.0, maxevals = 2
    );
    ST = typeof(tmp);
    discrete_fields = ST[
        discrete_system(ωs[i], psource, material;
            basis_order = basis_orders[i],
            basis_field_order = basis_field_orders[i],
            rtol = rtol, maxevals = maxevals
        )
    for i in eachindex(ωs)];

    # avoid right next to the surface due the boundary layer
    rs = 0.0:0.01:(R - 2*outer_radius(s1));
    xs = [ radial_to_cartesian_coordinates([r,0.2,1.2]) for r in rs];

    errors = [
        norm.(discrete_fields[i].coefficient_field.(xs) - pscat_fields[i].(xs)) ./ norm.(pscat_fields[i].(xs))
    for i in eachindex(ωs)];

    @test maximum(mean.(errors)) < 0.01
    @test maximum(maximum.(errors)) < 0.02

    # The scattering coefficients from the whole sphere has a smaller error. Indicating that the errors above occur more on the higher order modes which are smaller
    mat_coefs_pwaves = material_scattering_coefficients.(pwavemodes_azi);
    mat_coefs_discretes = material_scattering_coefficients.(discrete_fields;
        rtol = rtol,
        maxevals = maxevals
    );

    errors = [
        abs(norm(mat_coefs_pwaves[i][1:length(mat_coefs_discretes[i])]) / norm(mat_coefs_discretes[i]) - 1.0)
    for i in eachindex(ωs)];

    @test errors[1] < 1e-4
    @test errors[2] < 3e-3

## Radially symmetric scattering from a sphere

   # Define two radially symmetric sources, but only dispatch for radial symmetry on one of them
    sourceradial =  regular_spherical_source(medium, [1.0+0.0im];
       position = [0.0,0.0,0.0], symmetry = RadialSymmetry{3}()
    );

    sourceazi =  regular_spherical_source(medium, [1.0+0.0im];
       position = [0.0,0.0,0.0], symmetry = AzimuthalSymmetry{3}()
    );

    wavemodes_azi = [
        WaveMode(ωs[i], keffs[i], sourceazi, material;
           basis_order = basis_orders[i],
           basis_field_order = basis_field_orders[i]
        )
    for i in eachindex(ωs)];

    # Note there is no basis_field_order for RadialSymmetry below, due to symmetry restrictions
    wavemodes_radial = [
        WaveMode(ωs[i], keffs[i], sourceradial, material;
           basis_order = basis_orders[i]
        )
    for i in eachindex(ωs)];

    scat_fields_azi = scattering_field.(wavemodes_azi);
    scat_fields_radial = scattering_field.(wavemodes_radial);

    # pscat_field = scattering_field(pwavemode)

    # res = discrete_system_residue(pscat_field, ω, source, material, AzimuthalSymmetry{3}();
    #     basis_order = basis_order, mesh_points = 5,
    #     rtol = 1e-2, maxevals = Int(1e4)
    # )

    rtol = 1e-2; maxevals = Int(1e4);
    # this below is just to get the typeof tmp, to then avoid a weird unionall error which should hopefully go away when updating Julia at some point
    tmp = discrete_system(ωs[1], sourceradial, material;
        basis_order = 0, basis_field_order = 0, rtol = 1.0, maxevals = 4
    );
    ST = typeof(tmp);
    discrete_fields = ST[
        discrete_system(ωs[i], sourceradial, material;
            basis_order = basis_orders[i],
            basis_field_order = basis_field_orders[i],
            rtol = rtol, maxevals = maxevals
        )
    for i in eachindex(ωs)];

    xs = [ radial_to_cartesian_coordinates([r,0.0,0.0]) for r in rs];

    eff_scats_azi = [s.(xs) for s in scat_fields_azi];
    eff_scats_radial = [s.(xs) for s in scat_fields_radial];
    discrete_scats = [d.coefficient_field.(xs) for d in discrete_fields];

    # assuming radial or azimuthal symmetry should lead to exactly the same fields
    errors = [
        norm.(eff_scats_radial[i] - eff_scats_azi[i]) ./ norm.(eff_scats_radial[i])
    for i in eachindex(ωs)];

    @test maximum(maximum.(errors)) < 1e-12

    # the discrete method does have a difference
    errors = [
        norm.(eff_scats_radial[i] - discrete_scats[i]) ./ norm.(eff_scats_radial[i])
    for i in eachindex(ωs)];

    @test minimum(mean.(errors)) < 1e-4
    @test maximum(mean.(errors)) < 1e-3
    @test maximum(maximum.(errors)) < 2e-3

    mat_coefs_radial = material_scattering_coefficients.(wavemodes_radial);

    mat_coefs_disc_radial = material_scattering_coefficients.(discrete_fields;
        rtol = rtol,
        maxevals = maxevals
    );

    errors = norm.(mat_coefs_disc_radial - mat_coefs_radial) ./ norm.(mat_coefs_radial)

    @test errors[1] < 1e-5
    @test errors[2] < 2e-3

# Test the reduced radial discrete method

    pair_corr_inf(z) = hole_correction_pair_correlation([0.0,0.0,0.0],s1, [0.0,0.0,z],s1)

    polynomial_order = 15
    pair_corr_inf_smooth = smooth_pair_corr_distance(
        pair_corr_inf, a12;
        smoothing = 0.3, max_distance = 2R,
        polynomial_order = polynomial_order
    )

    # using Plots
    # zs = 0.0:0.01:(2R)
    # plot(pair_corr_inf_smooth,zs)
    # plot!(pair_corr_inf,zs, linestyle=:dash)

    gls_radial = gls_pair_radial_fun(pair_corr_inf_smooth, a12;
        sigma_approximation = false,
        polynomial_order = polynomial_order
    )

    discrete_field_radials = [
        discrete_system_radial(ωs[i], sourceradial, material, Symmetry(sourceradial,material);
            basis_order = basis_orders[i],
            basis_field_order = basis_field_orders[i],
            # polynomial_order = 15,
            # pair_corr_distance = pair_corr_inf
            gls_pair_radial = gls_radial
        )
    for i in eachindex(ωs)];

    discrete_rad_scats = [
        d.coefficient_field.(xs)
    for d in discrete_field_radials];

    errors = [
        norm.(discrete_rad_scats[i] - discrete_scats[i]) ./ norm.(discrete_scats[i])
    for i in eachindex(ωs)];

    @test mean(mean.(errors)) < 1e-3
    @test maximum(mean.(errors)) < 1e-3
    @test maximum(maximum.(errors)) < 4e-3

    mat_coefs_disc_radial2 = [
        material_scattering_coefficients(discrete_field_radials[i];
            rtol = rtol,
            maxevals = maxevals
        )
    for i in eachindex(ωs)];

    errors = norm.(mat_coefs_disc_radial2 - mat_coefs_disc_radial) ./ norm.(mat_coefs_disc_radial)
    @test maximum(errors) < 1e-3

# Calculate low frequency scattering
    # Linc = basis_orders[1] + basis_field_orders[1]
    # source_coefficients = regular_spherical_coefficients(sourceradial)(Linc,zeros(3),ωs[1]);
    #
    # material_low = Material(
    #     Sphere(outer_radius(material.shape) - outer_radius(s1)),
    #     species
    # );
    # #
    # effective_sphere = Particle(eff_medium, material_low.shape);
    # Tmat = MultipleScattering.t_matrix(effective_sphere, medium, ωs[1], Linc);
    # scat_coef_low = Tmat * source_coefficients;
end
