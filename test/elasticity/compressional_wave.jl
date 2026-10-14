
@testset "3D elastic compressional wave at direct incidence" begin

    symmetry = PlanarAzimuthalSymmetry{3}()

    # the secant method, to find a wavenumber where the eigensystem is singular
    function secant(f, x0)
        x1 = x0 * (1 + 1e-6); f0 = f(x0); f1 = f(x1)
        for _ = 1:80
            x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
            x0, f0 = x1, f1
            x1, f1 = x2, f(x2)
            if abs(x1 - x0) < 1e-14 * abs(x1) break end
        end
        return x1
    end

## Compare the eigensystem with the equations of the paper

    # Below is the part with azimuthal symmetry of equations (67a) and (67b) of [P. A. Martin and V. J. Pinfield, Wave Motion 134 (2025)], for one specie, written in the notation of that paper. For this symmetry F² = 0, and F³ is eliminated. The T-matrix used is not scaled.
    ω = 1.0
    medium = Elastic(3; ρ = 1.0, cp = complex(1 / 0.4), cs = 1 / (1.2 * (1 + im)))
    particle_medium = Elastic(3; ρ = 3.95, cp = complex(1 / 0.06), cs = complex(1 / 0.1))

    a = 1.0
    b = 2a
    n0 = 0.02
    N = 2

    kc = ω / medium.cp
    ks = ω / medium.cs

    # The block of the T-matrix for the pressure potential φ, and the shear potential Φ, of the degree l. There is one block of 3 potentials for each of the l^2 spherical harmonics with a degree lower than l.
    Tmat = t_matrix(Particle(particle_medium, Sphere(a)), medium, ω, N)
    Tblock(l) = Matrix(Tmat[3l^2 .+ (1:2), 3l^2 .+ (1:2)])

    γ(n) = sqrt((2n+1) / 4pi)
    # the Gaunt coefficient of equation (36) for the orders m = μ = 0
    gaunt(ν, n, q) = real(gaunt_coefficient(ν,0,n,0,q,0) * (1.0im)^(ν-n-q)) / 4pi
    # equation (62)
    Nq(q, x, y) = im * y * (x * diffsbesselj(q,x) * shankelh1(q,y) - y * sbesselj(q,x) * diffshankelh1(q,y))

    function paper_eigensystem(ξ; mode_conversion = true)
        # equation (67a) is F¹ = Ac * 𝔽¹, and (67b) is F³ = As * 𝔽³, where 𝔽 are the scattering coefficients (66), and the coefficient of (43b) is the Gaunt coefficient multiplied by (ν*(ν+1) + n*(n+1) - q*(q+1)) / (2n*(n+1))
        Ac = [
            n0 * 4pi * im / (kc * (kc^2 - ξ^2)) * 4pi * (1.0im)^(n-ν) *
            sum(γ(q) * Nq(q, ξ*b, kc*b) * gaunt(ν,n,q) for q = abs(n-ν):2:(n+ν))
        for n = 0:N, ν = 0:N]

        As = [
            (n == 0 || ν == 0) ? 0.0im :
            n0 * 4pi * im / (ks * (ks^2 - ξ^2)) * 4pi * (1.0im)^(n-ν) *
            sum(γ(q) * Nq(q, ξ*b, ks*b) * gaunt(ν,n,q) * (ν*(ν+1) + n*(n+1) - q*(q+1)) / (2n*(n+1)) for q = abs(n-ν):2:(n+ν))
        for n = 0:N, ν = 0:N]

        # The potentials used for the T-matrix, and the vector waves of the paper, differ by a factor which cancels in the products below
        T11 = Diagonal([Tblock(l)[1,1] for l = 0:N])
        T13 = Diagonal([Tblock(l)[1,2] for l = 0:N])
        T31 = Diagonal([Tblock(l)[2,1] for l = 0:N])
        T33 = Diagonal([Tblock(l)[2,2] for l = 0:N])

        # from (66) we have 𝔽¹ = T11 * F¹ + T13 * F³ and 𝔽³ = T31 * F¹ + T33 * F³, which together with F³ = As * 𝔽³ gives F³ in terms of F¹
        return mode_conversion ? I - (T11 + T13 * ((I - As * T33) \ (As * T31))) * Ac : I - T11 * Ac
    end

    specie = Specie(particle_medium, Sphere(a); volume_fraction = n0 * 4pi * a^3 / 3)
    micro = Microstructure(medium, [specie])

    MM = eigensystem(ω, micro, symmetry; basis_order = N)
    MM_pressure = eigensystem(ω, micro, symmetry; basis_order = N, mode_conversion = false)

    ξ = 0.47 + 0.11im
    @test maximum(abs.(MM(ξ) - paper_eigensystem(ξ))) < 1e-11
    @test maximum(abs.(MM_pressure(ξ) - paper_eigensystem(ξ; mode_conversion = false))) < 1e-11

    # the mode conversion is not negligible
    @test maximum(abs.(MM(ξ) - MM_pressure(ξ))) > 1e-2

## The eigensystem is singular at the effective wavenumber

    # The formula for the effective wavenumber is correct up to second order in the number density ε, so k_eff^2 differs by order ε^3 from where the eigensystem is singular. The first case is for two species, and the second for fine and coarse particles in a viscous liquid, where the shear waves need to be scaled.
    rs = [1.0, 0.6]
    species1(ε) = [
        Specie(particle_medium, Sphere(rs[i]); volume_fraction = ε * [1.0, 1.75][i] * 4pi * rs[i]^3 / 3)
    for i in eachindex(rs)]

    ρ = 1.0; c = 1.0;
    δ = 5e-3 * c / ω # the viscous skin depth
    ν = ω * δ^2 / 2 # the kinematic viscosity
    liquid = Elastic(3; ρ = ρ, cp = complex(c), cs = sqrt(-im * ω * ν))
    solid = Elastic(3; ρ = 3.95, cp = complex(6.7), cs = complex(3.7))

    species2(ε) = [
        Specie(solid, Sphere(2000 * δ); volume_fraction = 0.6 * ε),
        Specie(solid, Sphere(2 * δ); volume_fraction = 0.4 * ε)
    ]

    for (host, species) in [(medium, species1), (liquid, species2)]
        errors = map([1e-2, 5e-3]) do ε
            k_eff = wavenumber_compressional_low_volumefraction(ω, host, species(ε); basis_order = 2)
            MM = eigensystem(ω, Microstructure(host, species(ε)), symmetry; basis_order = 2)
            k_singular = secant(k -> det(MM(k)), k_eff)

            abs(k_singular^2 - k_eff^2)
        end

        @test 7 < errors[1] / errors[2] < 9
    end

    # there is one eigenvector, and no warning, though the eigensystem is not exactly singular for the wavenumber from the formula
    ε = 5e-2
    k_eff = wavenumber_compressional_low_volumefraction(ω, liquid, species2(ε))
    micro = Microstructure(liquid, species2(ε))

    eigvectors = @test_logs min_level=Base.CoreLogging.Warn eigenvectors(ω, k_eff, micro, symmetry; basis_order = 2)
    @test size(eigvectors) == (9, 2, 1)

    # the eigenvectors only have the spherical harmonics of the order m = 0
    ls, ms = spherical_harmonics_indices(2)
    @test norm(eigvectors[ms .!= 0, :, :]) == 0
    @test norm(eigensystem(ω, micro, symmetry; basis_order = 2)(k_eff) * vec(eigvectors[ms .== 0, :, 1])) / norm(eigvectors) < 10 * ε^2

## Reflection and transmission for fine particles in a viscous liquid

    # Particles which are much smaller than the viscous wavelength move with the liquid. Then the particulate material is the same as a homogeneous medium, where the density is the volume average of the densities, and the compressibility is the volume average of the compressibilities. This relies on the mode conversion: when it is neglected the relative errors of the reflection coefficients below are 0.12 and 0.07.
    δ = 2e-2 * c / ω # the viscous skin depth
    ν = ω * δ^2 / 2 # the kinematic viscosity
    liquid = Elastic(3; ρ = ρ, cp = complex(c), cs = sqrt(-im * ω * ν))

    solids = [
        Elastic(3; ρ = 4.0, cp = complex(7.0), cs = complex(4.0)),
        Elastic(3; ρ = 2.0, cp = complex(3.0), cs = complex(1.5))
    ]
    bulk_modulus(m) = m.ρ * (m.cp^2 - 4 * m.cs^2 / 3)

    # the pressure potential is treated as an acoustic wave
    psource = PlaneSource(Acoustic(3; ρ = ρ, c = c), [0.0, 0.0, 1.0])
    k = ω / c

    normal = [0.0, 0.0, -1.0] # an outward normal to both surfaces of the plate
    width = 4.0
    plate = Plate(normal, width, [0.0, 0.0, width / 2])
    halfspace = Halfspace(normal)

    for scale in [1.0, 0.5]
        volfracs = scale .* [0.06, 0.04]
        species = [
            Specie(solids[1], Sphere(2.0e-2 * δ); volume_fraction = volfracs[1]),
            Specie(solids[2], Sphere(1.4e-2 * δ); volume_fraction = volfracs[2])
        ]
        r = maximum(outer_radius.(species))

        ρ_eff = ρ * (1 - sum(volfracs)) + sum(volfracs .* [s.ρ for s in solids])
        β_eff = 1 / ((1 - sum(volfracs)) / (ρ * c^2) + sum(volfracs ./ bulk_modulus.(solids)))
        eff_medium = Acoustic(3; ρ = ρ_eff, c = sqrt(β_eff / ρ_eff))

        # The faces of the homogeneous plate are a distance r inside the faces of the particulate plate. The reflection and transmission coefficients of both plates are multiplied by exp(im * k * z), where z is the position of the first face for reflection, and of the second face for transmission, which leads to the phases below.
        plate_low = Plate(normal, width - 2r, [0.0, 0.0, width / 2])

        amps = planewave_amplitudes(ω, psource, eff_medium, plate_low)
        R_low = amps[1] * exp(-im * k * r)
        T_low = amps[2] * exp(im * k * r)

        k_eff = wavenumber_compressional_low_volumefraction(ω, liquid, species)

        material = Material(liquid, plate, species)
        wavemodes = @test_logs min_level=Base.CoreLogging.Warn WaveMode(ω, k_eff, psource, material)
        R, T = reflection_transmission_coefficients(wavemodes, psource, material)

        # these errors are mostly due to the effective wavenumber, which gives the phase change across the plate
        @test abs(R - R_low) / abs(R_low) < 5e-4
        @test abs(T - T_low) / abs(T_low) < 5e-4

        material = Material(liquid, halfspace, species)
        wavemode = @test_logs min_level=Base.CoreLogging.Warn WaveMode(ω, k_eff, psource, material)
        R = reflection_coefficient(wavemode, psource, material)

        # the reflection coefficient of a homogeneous halfspace whose face is at z = r
        q = eff_medium.ρ * eff_medium.c / (ρ * c)
        R_low = (q - 1) / (q + 1) * exp(2im * k * r)

        @test abs(R - R_low) / abs(R_low) < 3e-4
        @test abs(abs(R) / abs(R_low) - 1) < 2e-5
    end

## Cases which are not implemented

    species = [Specie(solids[1], Sphere(2.0e-2 * δ); volume_fraction = 0.05)]
    material = Material(liquid, plate, species)
    k_eff = wavenumber_compressional_low_volumefraction(ω, liquid, species)

    # a source which is not at direct incidence
    oblique_source = PlaneSource(Acoustic(3; ρ = ρ, c = c), [sin(0.3), 0.0, cos(0.3)])
    @test_throws ErrorException WaveMode(ω, k_eff, oblique_source, material)

    # a source in a medium which is different from the medium of the particles
    other_source = PlaneSource(Acoustic(3; ρ = 2ρ, c = c), [0.0, 0.0, 1.0])
    @test_throws ErrorException WaveMode(ω, k_eff, other_source, material)
end
