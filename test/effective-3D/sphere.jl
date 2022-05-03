using EffectiveWaves, ClassicalOrthogonalPolynomials
using LinearAlgebra, Statistics, Test

@testset "Effective sphere scattering" begin

## Set parameters
    # include("../discretisation/utils.jl")
    particle_medium = Acoustic(3; ρ=0.01, c=0.01);
    particle_medium = Acoustic(3; ρ=10.0, c=10.0);
    medium = Acoustic(3; ρ=1.0, c=1.0);

    R = 5.0
    r = 1.0

    separation_ratio = 1.02

    kas = [0.01,0.2]
    ks = kas ./ r

    vol_fraction = 0.12

    basis_orders = Int.(round.(4. .* kas)) .+ 1
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

    tol = 1e-7

## define sources and material
    eff_medium = effective_medium(medium, species)

    psource = PlaneSource(medium, [0.0,0.0,1.0]);
    source = plane_source(medium; direction = [0.0,0.0,1.0])

    sourceradial =  regular_spherical_source(medium, [1.0+0.0im];
       position = [0.0,0.0,0.0], symmetry = RadialSymmetry{3}()
    );

    sourceazi =  regular_spherical_source(medium, [1.0+0.0im];
       position = [0.0,0.0,0.0], symmetry = AzimuthalSymmetry{3}()
    );

    region_shape = Sphere([0.0,0.0,0.0], R)
    material = Material(Sphere(R),species);

    eff_medium = effective_medium(medium, species; numberofparticles = material.numberofparticles)
    ks_low = ωs ./ eff_medium.c


    keff_arr = [
        wavenumbers(ωs[i], medium, species;
            # num_wavenumbers = 4,
            basis_order = basis_orders[i],
            tol = tol,
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
    rs = 0.0:0.1:(R - 2.1 * outer_radius(s1));
    xs = [ radial_to_cartesian_coordinates([r,0.2,1.2]) for r in rs];

    errors = [
        norm.(discrete_fields[i].coefficient_field.(xs) - pscat_fields[i].(xs)) ./ norm.(pscat_fields[i].(xs))
    for i in eachindex(ωs)];

    @test maximum(mean.(errors)) < 0.01
    @test maximum(maximum.(errors)) < 0.02

    # The scattering coefficients from the whole sphere has a smaller error. Indicating that the errors above occur more on the higher order modes which are weaker
    mat_coefs_pwaves = material_scattering_coefficients.(pwavemodes_azi);
    mat_coefs_discretes = material_scattering_coefficients.(discrete_fields;
        rtol = rtol,
        maxevals = maxevals
    );

    errors = [abs(norm(mat_coefs_pwaves[i][1:length(mat_coefs_discretes[i])]) / norm(mat_coefs_discretes[i]) - 1.0) for i in eachindex(ωs)];
    @test errors[1] < 2e-4
    @test errors[2] < 3e-3

## Radially symmetric scattering from a sphere

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

    # source = sourceradial
    # basis_order = basis_orders[1]
    # ω = ωs[1]
    # basis_field_order = basis_field_orders[1]
    # basis_field_order = 4
    # T = Float64

    pair_corr_inf(z) = hole_correction_pair_correlation([0.0,0.0,0.0],s1, [0.0,0.0,z],s1)

    # gls_pair_radial(2.0,6.0)

    discrete_field_radials = [
        discrete_system_radial(ωs[i], sourceradial, material, Symmetry(sourceradial,material);
            basis_order = basis_orders[i],
            basis_field_order = basis_field_orders[i],
            polynomial_order = 20,
            pair_corr_distance = pair_corr_inf
        )
    for i in eachindex(ωs)];

    discrete_rad_scats = [
        d.coefficient_field.(xs)
    for d in discrete_field_radials];

    # scat_list = discrete_field_radial.coefficient_field.(xs);

    # the discrete method does have a difference
    # errors =
    # [
        # norm.(scat_list - discrete_scats[1]) ./ norm.(scat_list)
    # for i in eachindex(ωs)];

    errors = [
        norm.(discrete_rad_scats[i] - discrete_scats[i]) ./ norm.(discrete_scats[i])
    for i in eachindex(ωs)];

    minimum(mean.(errors[2]))
    mean(mean.(errors[2]))
    maximum(mean.(errors[2]))
    maximum(maximum.(errors[2]))

    mat_coefs_disc_radial2 = [
        material_scattering_coefficients(discrete_field_radials[i];
            rtol = rtol,
            maxevals = maxevals
        )
    for i in eachindex(ωs)];

    norm.(mat_coefs_disc_radial2 - mat_coefs_disc_radial) ./ norm.(mat_coefs_disc_radial)

    # abs.(mat_coefs_disc_radial2[1] - mat_coefs_disc_radial[1][1]) / abs.(mat_coefs_disc_radial[1][1])

    scat1 = [s[1] for s in discrete_rad_scats[2]];
    scat2 = [s[3] for s in discrete_rad_scats[2]];

    df0s = [d[1] for d in discrete_scats[2]];
    df1s = [d[3] for d in discrete_scats[2]];

    zs = [x[3] for x in xs];
    zs - norm.(xs)

    using Plots
    plot(zs,abs.(df0s))
    plot!(zs,abs.(scat1))

    plot(zs,abs.(df1s))
    plot!(zs,abs.(scat2))

    # scatter(abs.(discrete_field_radial[1,:]))
    # scatter(abs.(discrete_field_radial[2,:]))
    #
    # eff_scats_radial[1]

# Calculate low frequency scattering
    # Linc = basis_order + basis_field_order
    # source_coefficients = regular_spherical_coefficients(source)(Linc,zeros(3),ω);
    #
    # material_low = Material(
    #     Sphere(outer_radius(material.shape) - outer_radius(s1)),
    #     species
    # );
    #
    # effective_sphere = Particle(eff_medium, material_low.shape);
    # Tmat = MultipleScattering.t_matrix(effective_sphere, medium, ω, Linc);
    # scat_coef_low = Tmat * source_coefficients;
end

@testset "Discrete radial system for sphere scattering" begin

## Set parameters
    # include("../discretisation/utils.jl")
    particle_medium = Acoustic(3; ρ=0.1, c=0.1);
    particle_medium = Acoustic(3; ρ=10.01, c=10.01);
    medium = Acoustic(3; ρ=1.0, c=1.0);

    R = 5.0
    r = 1.0

    separation_ratio = 1.02

    kas = [0.01,0.2]
    kas = [0.002]
    kas = [0.1]
    kas = [0.3]
    ks = kas ./ r

    vol_fraction = 0.05
    vol_fraction = 0.15
    # vol_fraction = 0.04

    tol = 1e-7

    basis_orders = Int.(round.(4. .* kas)) .+ 1
    basis_orders = Int.(round.(0.0 .* kas)) .+ 1
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

    region_shape = Sphere([0.0,0.0,0.0], R)
    material = Material(Sphere(R),species);

    sourceradial = regular_spherical_source(medium, [1.0+0.0im];
       position = [0.0,0.0,0.0], symmetry = RadialSymmetry{3}()
    );

    t_matrices = get_t_matrices(sourceradial.medium, material.species, ωs[1], 3);
    t_diags = diag.(t_matrices)

    a12 = 2.0 * s1.exclusion_distance * outer_radius(s1)
    rs = 0.0:0.02:(R - a12);
    xs = [ radial_to_cartesian_coordinates([r,0.0,0.0]) for r in rs];

    eff_medium = effective_medium(medium, species; numberofparticles = material.numberofparticles)
    ks_low = ωs ./ eff_medium.c

    keff_arr = [
        wavenumbers(ωs[i], medium, species;
            num_wavenumbers = 4,
            # mesh_points = 6,
            # mesh_size = 1.0/5.0,
            basis_order = basis_orders[i],
            tol = 1e-5,
            numberofparticles = material.numberofparticles
        )
    for i in eachindex(ωs)]

    keffs = [ks[1] for ks in keff_arr]

    # Note there is no basis_field_order for RadialSymmetry below, due to symmetry restrictions
    wavemodes_radial = [
        WaveMode(ωs[i], keffs[i], sourceradial, material;
           basis_order = basis_orders[i]
        )
    for i in eachindex(ωs)];

    scat_fields_radial = scattering_field.(wavemodes_radial);
    eff_scats_radial = [s.(xs) for s in scat_fields_radial];


## Use a smooth pair- correlation
    # polynomial_order = 80
    polynomial_order = 15

    pair_corr_inf(z) = hole_correction_pair_correlation([0.0,0.0,0.0],s1, [0.0,0.0,z],s1)

    pair_corr_inf_smooth = smooth_pair_corr_distance(
        pair_corr_inf, a12;
        smoothing = 0.4, max_distance = 2R,
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

    # pair_radial_smooth = pair_radial_fun(pair_corr_inf_smooth, a12;
    #     sigma_approximation = false,
    #     polynomial_order = polynomial_order
    # )
    P = Legendre()

    function pair_radial_smooth2(r1,r2,u)
        Pus = P[u, 1:(polynomial_order + 1)] .* (2 .* (0:polynomial_order) .+ 1) ./ (4pi)

        return sum(Pus .* gls_radial(r1,r2))
    end

    plot!(zs,pair_radial_smooth2.(zs,0.0,0.0))
    plot!(zs,pair_radial_smooth2.(zs ./ 2,zs ./ 2,-1.0), linestyle=:dot)
    # costhetas = -1.0:0.1:1.0
    #
    # plot(costhetas, pair_radial.(4.01,0.0,costhetas))

    function pair_corr_smooth(x1,s1,x2,s2)
        if norm(x1) < 1e-10 || norm(x2) < 1e-10
            pair_radial_smooth2(norm(x1),norm(x2),0.0)
        else
            pair_radial_smooth2(norm(x1),norm(x2),dot(x1,x2) / (norm(x1)*norm(x2)))
        end
    end

    # pair_corr_inf_smooth2(z) = pair_corr_smooth([0.0,0.0,-z/2.0],s1,[0.0,0.0,z/2.0],s1)
    # plot!(pair_corr_inf_smooth2,zs, linestyle=:dash)

    gls_radial_simple(r1,r2) = [exp(-4r1 - 4r2)]

    pair_corr_simple(x1,s1,x2,s2) = gls_radial_simple(norm(x1),norm(x2))[1]

## Reduced radial method
    # import EffectiveWaves: discrete_system_radial
    discrete_field_radials = [
        EffectiveWaves.discrete_system_radial(ωs[i], sourceradial, material, Symmetry(sourceradial,material);
            basis_order = basis_orders[i],
            basis_field_order = basis_field_orders[i],
            # polynomial_order = 0,
            polynomial_order = polynomial_order,
            # mesh_points = 16,
            # mesh_points = 24,
            # pair_corr_distance = pair_corr_inf,
            # pair_corr_distance = pair_corr_inf_smooth,
            # gls_pair_radial = gls_radial_simple,
            gls_pair_radial = gls_radial,
            sigma_approximation = false
        )
    for i in eachindex(ωs)];

    discrete_rad_scats = [
        d.coefficient_field.(xs)
    for d in discrete_field_radials];

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
            rtol = 2e-3, maxevals = Int(1e5),
            # pair_corr = pair_corr_simple
            # pair_corr = pair_corr_smooth
        )
    for i in eachindex(ωs)];

    discrete_scats = [
        d.coefficient_field.(xs)
    for d in discrete_fields];

    mat_dcoefs = [
        material_scattering_coefficients(discrete_fields[i]; rtol = 1e-3,maxevals = Int(5e3))
    for i in eachindex(ωs)];

    mat_coefs_radial = material_scattering_coefficients.(wavemodes_radial);

   i = 1;

   abs(mat_dcoefs_radial[i][1])
   abs(mat_dcoefs[i][1])
   abs(mat_coefs_radial[i][1])

   abs(mat_dcoefs_radial[i][1] - mat_dcoefs[i][1]) / abs(mat_dcoefs[i][1])

   abs(mat_dcoefs[i][1] - mat_coefs_radial[i][1]) / abs(mat_coefs_radial[i][1])

   abs(mat_dcoefs_radial[i][1] - mat_coefs_radial[i][1]) / abs(mat_coefs_radial[i][1])

   j = 1;
   eff_rad0s = [s[j] for s in eff_scats_radial[i]];

   df_rad0s = [s[j] for s in discrete_rad_scats[i]];
   # df_rad1s = [s[3] for s in discrete_rad_scats[2]];
   df0s = [d[j] for d in discrete_scats[i]];
   # df1s = [d[3] for d in discrete_scats[2]];

   zs = [x[3] for x in xs];
   zs - norm.(xs);

   # using Plots
   # fun = abs;
   # plot(zs,fun.(df0s))
   # plot!(zs,fun.(df_rad0s), linestyle = :dash)
   # plot!(zs,fun.(eff_rad0s), linestyle = :dot)

   # plot(zs,fun.(df1s))
   # plot!(zs,fun.(df_rad1s))
end
