
@testset "low volume fraction 3D elastic compressional wavenumber" begin

## Compare with a reference solution for two species

    # The coefficients δ1 and δ2 of the number density, and of the number density squared, were calculated to 60 digits directly from equations (75), (94) and (100) of [P. A. Martin and V. J. Pinfield, Wave Motion 134 (2025)], with an independent T-matrix. They were also checked against the roots of the dispersion equation (67) of the same paper.
    ω = 1.0
    medium = Elastic(3; ρ = 1.0, cp = complex(1 / 0.4), cs = 1 / (1.2 * (1 + im)))
    particle_medium = Elastic(3; ρ = 3.95, cp = complex(1 / 0.06), cs = complex(1 / 0.1))

    δ1 = 1.0115357637519046 + 0.7638394897509469im
    δ2 = 0.5877715519008906 - 1.9319440852105250im

    ε = 1e-2
    rs = [1.0, 0.6]
    numdensities = ε .* [1.0, 1.75]
    species = [
        Specie(particle_medium, Sphere(rs[i]); volume_fraction = numdensities[i] * 4pi * rs[i]^3 / 3)
    for i in eachindex(rs)]

    k_eff = wavenumber_compressional_low_volumefraction(ω, medium, species; basis_order = 2)
    kp = ω / medium.cp

    @test k_eff^2 ≈ kp^2 + ε * δ1 + ε^2 * δ2 rtol = 1e-10

## Particles much smaller than the viscous wavelength move with the liquid

    # In this case the effective density is the volume average of the densities, and the effective compressibility is the volume average of the compressibilities, with no further terms. This relies on the mode conversion, which is 0.03 of k_eff^2 below.
    ρ = 1.0; c = 1.0;
    δ = 2e-2 * c / ω # the viscous skin depth
    ν = ω * δ^2 / 2 # the kinematic viscosity
    liquid = Elastic(3; ρ = ρ, cp = complex(c), cs = sqrt(-im * ω * ν))

    solids = [
        Elastic(3; ρ = 4.0, cp = complex(7.0), cs = complex(4.0)),
        Elastic(3; ρ = 2.0, cp = complex(3.0), cs = complex(1.5))
    ]
    volfracs = [0.06, 0.04]
    species = [
        Specie(solids[1], Sphere(2.0e-2 * δ); volume_fraction = volfracs[1]),
        Specie(solids[2], Sphere(1.4e-2 * δ); volume_fraction = volfracs[2])
    ]

    k_eff = wavenumber_compressional_low_volumefraction(ω, liquid, species)

    # the shear modulus of the liquid is negligible, so its bulk modulus is ρ * c^2
    bulk_modulus(m) = m.ρ * (m.cp^2 - 4 * m.cs^2 / 3)
    ρ_eff = ρ * (1 - sum(volfracs)) + sum(volfracs .* [s.ρ for s in solids])
    β_eff = 1 / ((1 - sum(volfracs)) / (ρ * c^2) + sum(volfracs ./ bulk_modulus.(solids)))

    @test abs(k_eff^2 - ω^2 * ρ_eff / β_eff) / abs(k_eff^2) < 5e-4

## Fine and coarse particles in a viscous liquid

    # The radius of the fine particles is 2 times the viscous skin depth, and the radius of the coarse particles is 2000 times. For the coarse particles the conversion of a pressure wave into a shear wave is of the order 1e867, which is too large for Float64, so this tests that the shear waves are scaled. The coefficients δ1 and δ2 of ε, and of ε squared, were calculated with 4096 bits of precision, and no scaling, from the same equations as for the first test, with the closed forms of the spherical Bessel functions. The same calculation gives δ1 and δ2 of the first test to 16 digits.
    δ = 5e-3 * c / ω # the viscous skin depth
    ν = ω * δ^2 / 2 # the kinematic viscosity
    liquid = Elastic(3; ρ = ρ, cp = complex(c), cs = sqrt(-im * ω * ν))
    solid = Elastic(3; ρ = 3.95, cp = complex(6.7), cs = complex(3.7))

    δ1 = 0.20797417440177765 + 0.17745622860119267im
    δ2 = 0.004186538281623636 + 0.032961020638257954im

    ε = 1e-2
    species = [
        Specie(solid, Sphere(2000 * δ); volume_fraction = 0.6 * ε),
        Specie(solid, Sphere(2 * δ); volume_fraction = 0.4 * ε)
    ]

    k_eff = wavenumber_compressional_low_volumefraction(ω, liquid, species; basis_order = 2)
    kp = ω / liquid.cp

    @test k_eff^2 ≈ kp^2 + ε * δ1 + ε^2 * δ2 rtol = 1e-10
end
