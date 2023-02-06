
"""
    scattering_field(args)

Returns a function which gives the average scattering coefficients for any vector `x` inside the material. This field is defined by Equation (3.13) in [AL Gower and G Kristensson, "Effective waves for random three-dimensional particulate materials", (2021)](https://arxiv.org/pdf/2010.00934.pdf)
"""
scattering_field


"Calculates the effective wavenumbers and return Vector{EffectivePlaneWaveMode}."
function WaveModes(ω::T, source::AbstractSource, material::Material{S}; kws...) where {T,S<:Shape} # without the parametric types we get a "Unreachable reached" error

    # The wavenumbers are calculated without knowledge of the material symmetry. This is because the plane-wave symmetry leads to all possible wavenumbers and is simple to calculate.
    k_effs = wavenumbers(ω, material.microstructure; kws... )

    # The wavemodes need to know the material symmetry as the eigenvectors do depend on material shape and symetry.
    wave_effs = [
        WaveMode(ω, k_eff, source, material; kws...)
    for k_eff in k_effs]

    return wave_effs
end

"""
    WaveMode(ω::T, wavenumber::Complex{T}, eigenvectors::Array{Complex{T}}, ::SetupSymmetry; kws...)

Returns a concrete subtype of AbstractWaveMode depending on the SetupSymmetry. The returned type should have all the necessary fields to calculate scattered waves (currently not true for EffectivePlanarWaves).
"""
function WaveMode(ω::T, wavenumber::Complex{T}, source::AbstractSource, material::Material; kws...) where T

    eigvectors = eigenvectors(ω, wavenumber, source, material; kws...)

    α = solve_boundary_condition(ω, wavenumber, eigvectors, source, material; kws...)

    # After this normalisation, sum(eigvectors, dims = 3) will satisfy the boundary conditions
    eigvectors = [eigvectors[i] * α[i[3]] for i in CartesianIndices(eigvectors)]

    return EffectiveRegularWaveMode(ω, wavenumber, source, material, eigvectors; kws...)
end

function WaveMode(ω::T, wavenumber::Complex{T}, psource::PlaneSource{T,Dim,1}, material::Material{Halfspace{T,Dim}};
    tol::T = 1e-6, kws...) where {T,Dim}

    direction = transmission_direction(wavenumber, (ω / psource.medium.c) * psource.direction, material.shape.normal)
    eigvectors = eigenvectors(ω, wavenumber, psource, material; direction_eff = direction, kws...)

    α = solve_boundary_condition(ω, wavenumber, eigvectors, psource, material; kws...)

    # After this normalisation, sum(eigvectors, dims = 3) will satisfy the boundary conditions
    eigvectors = [eigvectors[i] * α[i[3]] for i in CartesianIndices(eigvectors)]

    return EffectivePlaneWaveMode(ω, wavenumber, direction, eigvectors)
end

function WaveMode(ω::T, wavenumber::Complex{T}, psource::PlaneSource{T,Dim,1}, material::Material{Plate{T,Dim}}; kws...) where {T,Dim}

    if psource.medium != material.microstructure.medium @error mismatched_medium end

    # First we calculate the outward pointing normal
    n = material.shape.normal;
    n = - n .* sign(real(dot(n,psource.direction)));

    k = ω / psource.medium.c

    direction1 = transmission_direction(wavenumber, k * psource.direction, n)
    eigvectors1 = eigenvectors(ω, wavenumber, psource, material; direction_eff = direction1, kws...)

    # we choose direction2 so that k2 .* direction2 = - k1 .* direction1, where k1 = wavenumber, and k2 = - wavenumber
    direction2 = direction1
    eigvectors2 = eigenvectors(ω, - wavenumber, psource, material; direction_eff = direction2, kws...)

    α = solve_boundary_condition(ω, wavenumber, eigvectors1, eigvectors2, psource, material; kws...)

    # apply normalisation
    eigvectors1 = eigvectors1 .* α[1]
    eigvectors2 = eigvectors2 .* α[2]

    mode1 = EffectivePlaneWaveMode(ω, wavenumber, direction1, eigvectors1)
    mode2 = EffectivePlaneWaveMode(ω, - wavenumber, direction2, eigvectors2)

    return [mode1,mode2]
end

function solve_boundary_condition(ω::T, wavenumber::Complex{T}, eigvectors::Array, source::AbstractSource, material::Material; kws...) where T
    return solve_boundary_condition(ω, wavenumber, eigvectors, source, material, Symmetry(source,material); kws...)
end

function solve_boundary_condition(ω::T, wavenumber::Complex{T}, eigvectors1::Array, eigvectors2::Array, source::AbstractSource, material::Material; kws...) where T
    return solve_boundary_condition(ω, wavenumber, eigvectors1, eigvectors2, source, material, Symmetry(source,material); kws...)
end

"""
    convert_eigenvector_basis(medium::PhysicalMedium,sym::AbstractSymmetry,eigvecs::Array)

The eigenvectors from high symmetric scenarios are smaller then the more general scenarios. This function pads with zeros the more symmetric cases to match more general cases, so that we can use the functions for both.
"""
convert_eigenvector_basis(medium::PhysicalMedium,sym::AbstractSymmetry,eigvecs::Array) = eigvecs

function eigenvectors(ω::T, k_eff::Complex{T}, source::AbstractSource, material::Material; kws...) where T<:AbstractFloat
    eigenvectors(ω, k_eff::Complex{T}, material.microstructure, Symmetry(source,material);
        kws...
    )
end

# For plane waves, it is simpler to write all cases in the format for the most general case. For example, for PlanarAzimuthalSymmetry the eignvectors are much smaller. So we will turn these into the more general eigvector case by padding it with zeros.
# function eigenvectors(ω::T, k_eff::Complex{T}, source::PlaneSource{T}, material::Material{S}; kws...) where {T<:AbstractFloat,Dim,S<:Union{Plate,Halfspace}}
#
#     eigvecs = eigenvectors(ω, k_eff, source.medium, material.microstructure.species, Symmetry(source,material); kws...)
#
#     if Symmetry(source,material) == PlanarAzimuthalSymmetry{Dim}()
#         eigvecs = azimuthal_to_planar_eigenvector(typeof(source.medium),eigvecs)
#     end
#
#     return eigvecs
#
# end

function eigenvectors(ω::T, k_eff::Complex{T}, micro::Microstructure, symmetry::AbstractSymmetry;
        tol::T = 1e-4, kws...
    ) where T<:AbstractFloat

    MM = eigensystem(ω, micro, symmetry; kws...)

    # calculate eigenvectors
    MM_svd = svd(MM(k_eff))
    inds = findall(MM_svd.S .< tol)

    if isempty(inds)
        @warn("No eigenvectors found with the tolerance tol = $tol. Will return only one eigenvector with the tolernace $(MM_svd.S[end])")
        inds = [length(MM_svd.S)]
    end

    #NOTE: MM(k_eff) ≈ MM_svd.U * diagm(0 => MM_svd.S) * MM_svd.Vt
    eigvectors = MM_svd.V[:,inds]

    # Reshape to separate different species and eigenvectors
    S = length(micro.species)
    eigvectors = reshape(eigvectors,(:,S,size(eigvectors,2)))

    # pads with zeros if necessary to match the more general case with less symmetry
    eigvectors = convert_eigenvector_basis(micro.medium,symmetry,eigvectors)

    return eigvectors
end
