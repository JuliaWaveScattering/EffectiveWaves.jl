# The effective compressional wave in an elastic medium filled with spherical particles, for a plane compressional wave at direct incidence on a plate or halfspace.
#
# At direct incidence everything is symmetric about the direction of the incident wave, so only the spherical waves with m = 0 are present. In this case the scattered shear waves, when averaged over the positions of the particles, do not form a plane wave. So the boundary conditions, and the reflected and transmitted waves, have the same form as for acoustics when written in terms of the pressure potential φ, where the displacement of the compressional wave is ∇φ. The shear waves only appear in the eigensystem, where they change how a particle scatters a pressure wave into a pressure wave.

import ElasticWaves: t_matrix_shear_scaled

"""
    eigensystem(ω, micro::ParticulateMicrostructure{3,Elastic{3,T}}, ::PlanarAzimuthalSymmetry{3}; basis_order = 2, mode_conversion = true)

Returns a function `MM` where `MM(k_eff)` is the matrix of the eigensystem for the effective compressional wave in an elastic medium filled with spherical particles, for a plane compressional wave at direct incidence. When `k_eff` is the effective wavenumber, such as given by [`wavenumber_compressional_low_volumefraction`](@ref), the null space of `MM(k_eff)` gives the average scattering coefficients of the pressure potential φ, where the displacement of the compressional wave is `∇φ`. These have the same form as the eigenvectors for an `Acoustic` medium, so they can be used with the methods which calculate reflection and transmission for acoustics.

The eigensystem is the part with azimuthal symmetry of equations (67a) and (67b) in [P. A. Martin and V. J. Pinfield, "Elastodynamic multiple scattering: effective wavenumbers in three-dimensional elastic media", Wave Motion 134 (2025)], after eliminating the shear waves, and written in terms of the scattering coefficients given by equation (66). The result is the eigensystem for acoustics with the T-matrix

    Tpp + Tps * inv(I - As(k_eff) * Tss) * As(k_eff) * Tsp

where `Tpp` scatters a pressure wave into a pressure wave, `Tsp` converts a pressure wave into a shear wave, `As(k_eff)` carries a shear wave from one particle to another, `Tss` scatters a shear wave into a shear wave, and `Tps` converts a shear wave back into a pressure wave. The inverse accounts for the shear wave being scattered any number of times between the particles.

Use `mode_conversion = false` to neglect the shear waves between the particles, which is the eigensystem for acoustics with the T-matrix `Tpp`. The transmission through a plate changes little, but the reflection coefficient then has a relative error which is first order in the particle volume fraction.

The particles are assumed to not overlap, but to be otherwise uncorrelated, which is called the hole correction.

# Example

To calculate the reflection and transmission from a plate, use a plane wave source in the `Acoustic` medium which has the same density and pressure wave speed as the elastic medium:

```julia
medium = Elastic(3; ρ = 1.0, cp = 1.0, cs = 0.01 - 0.01im)
species = [Specie(Elastic(3; ρ = 4.0, cp = 7.0, cs = 4.0), Sphere(0.02); volume_fraction = 0.1)]

ω = 1.0
k_eff = wavenumber_compressional_low_volumefraction(ω, medium, species)

material = Material(medium, Plate([0.0,0.0,-1.0], 4.0, [0.0,0.0,2.0]), species)
psource = PlaneSource(Acoustic(3; ρ = medium.ρ, c = medium.cp), [0.0,0.0,1.0])

wavemodes = WaveMode(ω, k_eff, psource, material)
R, T = reflection_transmission_coefficients(wavemodes, psource, material)
```

Here `R` and `T` are the reflection and transmission coefficients of the pressure potential φ, which are the same as those of the normal stress. For the displacement the reflection coefficient is `-R`.
"""
function eigensystem(ω::T, micro::ParticulateMicrostructure{3,Elastic{3,T}}, ::PlanarAzimuthalSymmetry{3};
        basis_order::Int = 2,
        mode_conversion::Bool = true,
        kws...) where {T<:AbstractFloat}

    medium = micro.medium
    species = micro.species

    if length(micro.paircorrelations[1].r) > 1
        error("For an elastic medium the eigensystem is only implemented for the hole correction pair correlation.")
    end

    kp = ω / medium.cp
    ks = ω / medium.cs
    S = length(species)
    ho = basis_order
    len = (ho+1) * S

    # For each basis order l, Tpp scatters the pressure potential φ into φ, Tsp converts φ into the shear potential Φ, Tps converts Φ into φ, and Tss scatters Φ into Φ. Here Tsp and Tps are multiplied by exp(im * ks * a), and Tss is multiplied by exp(2im * ks * a), where a is the radius of the particle. Without these factors they are exponentially large for a shear wave which decays quickly, as it does in a viscous liquid. See t_matrix_pressure_scaled and t_matrix_shear_scaled.
    Ts = [t_matrix_pressure_scaled(s.particle, medium, ω, ho) for s in species]
    Tpp = vcat(first.(Ts)...)
    Tsp = vcat(last.(Ts)...)

    Ts = [t_matrix_shear_scaled(s.particle, medium, ω, ho) for s in species]
    Tps = vcat(first.(Ts)...)
    Tss = vcat(last.(Ts)...)

    rs = outer_radius.(species)
    as = [
        micro.paircorrelations[i,j].minimal_distance
    for i in eachindex(species), j in eachindex(species)]

    # The coefficients which translate a wave with basis order dl, centred at one particle, into a wave with basis order l, centred at another particle. These are zero when l + dl + l1 is odd.
    gaunts = [
        (abs(dl-l) <= l1 <= dl+l) ? sqrt(4 * T(pi) * (2l1+1)) * Complex{T}(im)^(-l1) * gaunt_coefficient(T,dl,0,l,0,l1,0) : zero(Complex{T})
    for l1 = 0:2ho, l = 0:ho, dl = 0:ho]

    function MM(keff::Complex{T})::Matrix{Complex{T}}

        # carries a pressure wave between two particles whose centres are at least a distance as[s1,s2] apart
        Np = [
            as[s1,s2] * kernelN3D(l1, kp*as[s1,s2], keff*as[s1,s2])
        for l1 = 0:2ho, s1 = 1:S, s2 = 1:S] ./ (keff^2 - kp^2)

        # the eigensystem for acoustics is I + Diagonal(T-matrix) * Kp
        Kp = [
            number_density(species[s2]) * sum(gaunts[l1+1,l+1,dl+1] * Np[l1+1,s1,s2] for l1 = abs(dl-l):2:(dl+l))
        for l = 0:ho, s1 = 1:S, dl = 0:ho, s2 = 1:S]
        Kp = reshape(Kp,len,len)

        if !mode_conversion
            return I + Diagonal(Tpp) * Kp
        end

        # Carries a shear wave between two particles. The same as Np with kp replaced by ks, except it is multiplied by exp(-im * ks * (rs[s1] + rs[s2])), which cancels the factors that multiply Tsp, Tps, and Tss. Without this the shear wave is exponentially small when it decays quickly. The exponential below equals 1 when the particles can touch.
        Ns = [
            - as[s1,s2] * kernel_shear_scaled(l1, keff*as[s1,s2], ks*as[s1,s2]) *
            exp(im * ks * (as[s1,s2] - rs[s1] - rs[s2]))
        for l1 = 0:2ho, s1 = 1:S, s2 = 1:S] ./ (keff^2 - ks^2)

        # the shear wave scattered from the particle s2, with basis order dl, which arrives at the particle s1 with basis order l. There are no shear waves for the basis order 0.
        As = [
            (l == 0 || dl == 0) ? zero(Complex{T}) :
            - number_density(species[s2]) / (2l * (l+1)) * sum(
                gaunts[l1+1,l+1,dl+1] * (l*(l+1) + dl*(dl+1) - l1*(l1+1)) * Ns[l1+1,s1,s2]
            for l1 = abs(dl-l):2:(dl+l))
        for l = 0:ho, s1 = 1:S, dl = 0:ho, s2 = 1:S]
        As = reshape(As,len,len)

        # A pressure wave is converted into a shear wave by Tsp. This shear wave is then scattered any number of times between the particles, which gives the inverse below, before it is converted back into a pressure wave by Tps.
        Tpsp = Diagonal(Tps) * ((I - As * Diagonal(Tss)) \ (As * Diagonal(Tsp)))

        return I + (Diagonal(Tpp) + Tpsp) * Kp
    end

    return MM
end

eigensystem(::T, ::ParticulateMicrostructure{3,Elastic{3,T}}, ::AbstractPlanarSymmetry{3}; kws...) where {T<:AbstractFloat} = error("For an elastic medium the eigensystem is only implemented for a compressional plane wave at direct incidence, which has PlanarAzimuthalSymmetry.")

# The effective wavenumber from wavenumber_compressional_low_volumefraction is correct up to second order in the particle volume fraction. Due to this, the smallest singular value of the eigensystem is of the order of the volume fraction squared, instead of zero, so a smaller tolerance can not be met. The factor of 10 below is because the smallest singular value is as large as 7 times the volume fraction squared for some particles. The tolerance only decides when to warn that k_eff is far from an effective wavenumber, as the eigenvector returned is always the one with the smallest singular value.
function eigenvectors(ω::T, k_eff::Complex{T}, micro::ParticulateMicrostructure{3,Elastic{3,T}}, symmetry::PlanarAzimuthalSymmetry{3};
        tol::T = 1e-4, kws...
    ) where T<:AbstractFloat

    tol = max(tol, 10 * volume_fraction(micro.species)^2)

    eigvectors = invoke(eigenvectors, Tuple{T,Complex{T},Microstructure,AbstractSymmetry}, ω, k_eff, micro, symmetry; tol = tol, kws...)

    # There is only one compressional wave for each effective wavenumber, which is the eigenvector with the smallest singular value. This is the last one.
    return eigvectors[:,:,end:end]
end

# The acoustic medium which has the same pressure wave as the elastic medium. A field in this acoustic medium represents the pressure potential φ.
pressure_medium(medium::Elastic{Dim,T}) where {Dim,T} = Acoustic(Dim; ρ = medium.ρ, c = medium.cp)

# Replaces the elastic medium of the material with pressure_medium, so that the methods for acoustics can be used.
function pressure_material(psource::PlaneSource{T,3,1,Acoustic{T,3}}, material::Material{Sh,<:ParticulateMicrostructure{3,Elastic{3,T}}}) where {T,Sh<:Shape}

    micro = material.microstructure
    medium = pressure_medium(micro.medium)

    if psource.medium != medium
        error("For an elastic medium the source needs to be in the medium Acoustic(3; ρ = medium.ρ, c = medium.cp), where medium is the elastic medium of the microstructure. " * mismatched_medium)
    end

    return Material(material.shape, ParticulateMicrostructure(medium, micro.species, micro.paircorrelations))
end

convert_eigenvector_basis(medium::Elastic{3,T}, sym::PlanarAzimuthalSymmetry{3}, eigvecs::Array{Complex{T}}) where T = convert_eigenvector_basis(pressure_medium(medium), sym, eigvecs)

## The eigenvectors are the average scattering coefficients of the pressure potential, and have the same form as for acoustics, so below we use the methods for acoustics.

function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, psource::PlaneSource{T,3,1,Acoustic{T,3}}, material::Material{Halfspace{T,3},<:ParticulateMicrostructure{3,Elastic{3,T}}}, sym::AbstractPlanarSymmetry{3};
        kws...
    ) where T

    return solve_boundary_condition(ω, k_eff, eigvectors, psource, pressure_material(psource, material), sym; kws...)
end

function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors1::Array{Complex{T}}, eigvectors2::Array{Complex{T}}, psource::PlaneSource{T,3,1,Acoustic{T,3}}, material::Material{Plate{T,3},<:ParticulateMicrostructure{3,Elastic{3,T}}}, sym::AbstractPlanarSymmetry{3};
        kws...
    ) where T

    return solve_boundary_condition(ω, k_eff, eigvectors1, eigvectors2, psource, pressure_material(psource, material), sym; kws...)
end

function reflection_coefficient(wavemode::EffectivePlaneWaveMode{T,Dim}, psource::PlaneSource{T,3,1,Acoustic{T,3}}, material::Material{Halfspace{T,3},<:ParticulateMicrostructure{3,Elastic{3,T}}}) where {T<:AbstractFloat,Dim}

    return reflection_coefficient(wavemode, psource, pressure_material(psource, material))
end

function reflection_transmission_coefficients(wavemodes::Vector{E}, psource::PlaneSource{T,3,1,Acoustic{T,3}}, material::Material{Plate{T,3},<:ParticulateMicrostructure{3,Elastic{3,T}}}) where {T<:AbstractFloat,Dim, E<:EffectivePlaneWaveMode{T,Dim}}

    return reflection_transmission_coefficients(wavemodes, psource, pressure_material(psource, material))
end
