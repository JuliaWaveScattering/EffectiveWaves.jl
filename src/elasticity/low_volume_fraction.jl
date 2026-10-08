export wavenumber_compressional_low_volumefraction

import SpecialFunctions: hankelh1x
import ElasticWaves: t_matrix_pressure_scaled

wavenumber_compressional_low_volumefraction(ω::Number, medium::Elastic{3}, specie::Specie; kws...) = wavenumber_compressional_low_volumefraction(ω, medium, [specie]; kws...)

wavenumber_compressional_low_volumefraction(ωs::AbstractVector{<:Number}, medium::Elastic{3}, species::Species; kws...) = [wavenumber_compressional_low_volumefraction(ω, medium, species; kws...) for ω in ωs]

wavenumber_compressional_low_volumefraction(ωs::AbstractVector{<:Number}, medium::Elastic{3}, specie::Specie; kws...) = wavenumber_compressional_low_volumefraction(ωs, medium, [specie]; kws...)

# The kernel which carries a shear wave between two particles whose centres are a distance b apart: y jₗ(x) hₗ₊₁(y) - x jₗ₊₁(x) hₗ(y), where x = kp * b and y = ks * b, multiplied by exp(-im*y). The scaling avoids underflow when the shear wave decays quickly, as it does in a viscous liquid.
function kernel_shear_scaled(l::Int, x::Complex{T}, y::Complex{T}) where T<:AbstractFloat
    scaled_shankelh1(n) = sqrt(pi / (2y)) * hankelh1x(n + T(1)/2, y)

    return y * sbesselj(l, x) * scaled_shankelh1(l+1) - x * sbesselj(l+1, x) * scaled_shankelh1(l)
end

"""
    wavenumber_compressional_low_volumefraction(ω::T, medium::Elastic{3,T}, species::Species{3}; basis_order::Int = 2)

Explicit formula for the effective compressional wavenumber `k_eff` of an elastic `medium` filled with spherical particles, based on a low particle volume fraction expansion:

    k_eff^2 = kp^2 + K1 + K2pp + K2sp,

where `kp = ω / medium.cp`, the term `K1` is first order in the number density of the particles, while `K2pp` and `K2sp` are second order. The formula is given by equations (75) and (100) in [P. A. Martin and V. J. Pinfield, "Elastodynamic multiple scattering: effective wavenumbers in three-dimensional elastic media", Wave Motion 134 (2025)], here written for many `species` of spherical particles.

The terms `K1` and `K2pp` have the same form as for acoustics, see [`wavenumber_low_volumefraction`](@ref), except they use the pressure to pressure part of the elastic T-matrix. The term `K2sp` is due to mode conversion: the pressure wave is converted into a shear wave at one particle, and then back into a pressure wave at another particle.

The particles are assumed to not overlap, but to be otherwise uncorrelated, which is called the hole correction.

A viscous liquid, with kinematic viscosity `ν`, is the same as an elastic medium with the complex shear wave speed `cs = sqrt(-im * ω * ν)`.
"""
function wavenumber_compressional_low_volumefraction(ω::T, medium::Elastic{3,T}, species::Species{3};
        basis_order::Int = 2, verbose::Bool = true
    ) where T<:AbstractFloat

    volfrac = volume_fraction(species)
    if volfrac >= 0.4 && verbose
        @warn("the volume fraction $(round(100*volfrac))% is too high, expect a relative error of approximately $(round(100*volfrac^3.0))%")
    end

    # background wavenumbers
    kp = ω / medium.cp
    ks = ω / medium.cs

    # For each basis order l the T-matrix has a 3 x 3 block which acts on the coefficients of the potentials [φ, Φ, χ]. We need only how the pressure potential φ scatters into φ, called Tpp, and how φ scatters into the shear potential Φ, called Tsp.
    # Below Tsp is multiplied by exp(im * ks * a), where a is the radius of the particle. Without this factor Tsp is exponentially large for a shear wave which decays quickly, as in a viscous liquid, and the T-matrix itself can not be calculated.
    Ts = [t_matrix_pressure_scaled(s.particle, medium, ω, basis_order) for s in species]

    Tpp = first.(Ts)
    Tsp = last.(Ts)

    numdensities = number_density.(species)
    rs = outer_radius.(species)

    # The minimal distance between the centres of two particles
    bs = [s1.separation_ratio * outer_radius(s1) + s2.separation_ratio * outer_radius(s2) for s1 in species, s2 in species]

    # W(l,dl,l1) = (2l+1) * (2dl+1) * (2l1+1) * wigner3j(l,dl,l1,0,0,0)^2, which is zero when l + dl + l1 is odd
    W(l::Int, dl::Int, l1::Int) = sqrt((2l+1) * (2dl+1) * (2l1+1) / (4 * T(pi))) * real(Complex{T}(im)^(l-dl-l1) * gaunt_coefficient(T,l,0,dl,0,l1,0))

    K1 = - im * 4pi / kp * sum(
        numdensities[s] * (2l + 1) * Tpp[s][l+1]
    for l = 0:basis_order, s in eachindex(species))

    K2pp = im * 8pi^2 / kp^3 * sum(
        W(l,dl,l1) * sum(
            bs[s1,s2] * d3D(kp * bs[s1,s2], l1) *
            Tpp[s1][l+1] * Tpp[s2][dl+1] * numdensities[s1] * numdensities[s2]
        for s1 in eachindex(species), s2 in eachindex(species))
    for l = 0:basis_order for dl = 0:basis_order for l1 = abs(l-dl):2:(l+dl))

    # The conversion from Φ back to φ equals l * (l+1) * (kp / ks) multiplied by the conversion from φ to Φ, because the T-matrix is reciprocal, so only Tsp appears below. It is not taken from the T-matrix itself, because for a shear wave which decays quickly, as in a viscous liquid, the part of the T-matrix for an incident shear wave overflows much sooner than Tsp does.
    # In the same case the shear wave between two particles, whose centres are a distance b apart, is exponentially small. The kernel below is multiplied by exp(-im * ks * b), and each Tsp by exp(im * ks * a), so that their product is evaluated without overflow. The exponential which is left equals 1 when the particles can touch.
    K2sp = if basis_order == 0
        zero(Complex{T})
    else
        im * 8pi^2 / (ks * (ks^2 - kp^2)) * sum(
            W(l,dl,l1) * (l*(l+1) + dl*(dl+1) - l1*(l1+1)) * sum(
                bs[s1,s2] * kernel_shear_scaled(l1, kp * bs[s1,s2], ks * bs[s1,s2]) *
                exp(im * ks * (bs[s1,s2] - rs[s1] - rs[s2])) *
                Tsp[s1][l+1] * Tsp[s2][dl+1] *
                numdensities[s1] * numdensities[s2]
            for s1 in eachindex(species), s2 in eachindex(species))
        for l = 1:basis_order for dl = 1:basis_order for l1 = abs(l-dl):2:(l+dl))
    end

    # effective wavenumber squared up to second order in the particle volume fraction
    kT2::Complex{T} = kp^2 + K1 + K2pp + K2sp

    kT = if abs(imag(sqrt(kT2)) / real(sqrt(kT2))) < eps(T) 
        real(sqrt(kT2))  < 0 ? -sqrt(kT2) : sqrt(kT2)
    else 
        imag(sqrt(kT2)) > 0 ? sqrt(kT2) : -sqrt(kT2)
    end

    return kT
end
