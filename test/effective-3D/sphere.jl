using EffectiveWaves
using LinearAlgebra, Statistics, Test


@testset "Effective sphere scattering" begin

    # include("test_discrete_solution.jl")
# Set parameters
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

# define sources and material
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


# [mat_coefs_pwave[inds] mat_coefs_disc2[inds] mat_coefs_disc[inds]]

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
