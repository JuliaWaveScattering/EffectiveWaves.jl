"Calculates the effective wavenumbers and return Vector{EffectivePlaneWaveMode}."
function WaveModes(ω::T, source::AbstractSource, material::Material{Dim,S,Sps}; kws...) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}} # without the parametric types we get a "Unreachable reached" error

    # The wavenumbers are calculated without knowledge of the materail symmetry. This is because the plane-wave symmetry leads to all possible wavenumbers and is simple to calculate.
    k_effs = wavenumbers(ω, source.medium, material.species; kws... )

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
function WaveMode(ω::T, wavenumber::Complex{T}, source::AbstractSource{T}, material::Material{Dim}; kws...) where {T,Dim}

    eigvectors = eigenvectors(ω, wavenumber, source, material; kws...)

    α = solve_boundary_condition(ω, wavenumber, eigvectors, source, material; kws...)

    # After this normalisation, sum(eigvectors, dims = 3) will satisfy the boundary conditions
    eigvectors = [eigvectors[i] * α[i[3]] for i in CartesianIndices(eigvectors)]

    return EffectiveRegularWaveMode(ω, wavenumber, source, material, eigvectors; kws...)
end

function WaveMode(ω::T, wavenumber::Complex{T}, psource::PlaneSource{T,Dim,1}, material::Material{Dim,Halfspace{T,Dim}};
    tol::T = 1e-6, kws...) where {T,Dim}

    direction = transmission_direction(wavenumber, (ω / psource.medium.c) * psource.direction, material.shape.normal)
    eigvectors = eigenvectors(ω, wavenumber, psource, material; direction_eff = direction, kws...)

    α = solve_boundary_condition(ω, wavenumber, eigvectors, psource, material; kws...)

    # After this normalisation, sum(eigvectors, dims = 3) will satisfy the boundary conditions
    eigvectors = [eigvectors[i] * α[i[3]] for i in CartesianIndices(eigvectors)]

    return EffectivePlaneWaveMode(ω, wavenumber, direction, eigvectors)
end

function WaveMode(ω::T, wavenumber::Complex{T}, psource::PlaneSource{T,Dim,1}, material::Material{Dim,Plate{T,Dim}};
    tol::T = 1e-6, kws...) where {T,Dim}

    direction = transmission_direction(wavenumber, (ω / psource.medium.c) * psource.direction, material.shape.normal)
    eigvecs1 = eigenvectors(ω, wavenumber, psource, material; direction_eff = direction, kws...)
    eigvecs2 = eigenvectors(ω, - wavenumber, psource, material; direction_eff = direction, kws...)

    α = solve_boundary_condition(ω, k_eff, eigvecs1, eigvecs2, psource, material; kws...)
    @error "not yet implemented. Needs two wave modes."
    # After this normalisation, sum(eigvectors, dims = 3) will satisfy the boundary conditions
    # eigvectors = [eigvectors[i] * α[i[3]] for i in CartesianIndices(eigvectors)]

    return EffectivePlaneWaveMode(ω, wavenumber, direction, eigvectors)
end

eigensystem(ω::T, source::AbstractSource{T}, material::Material; kws...) where T<:AbstractFloat = eigensystem(ω, source.medium, material.species, setupsymmetry(source,material); numberofparticles = material.numberofparticles, kws...)

function eigenvectors(ω::T, k_eff::Complex{T}, source::AbstractSource, material::Material;
        tol::T = 1e-4,
        kws...) where T<:Number

    MM = eigensystem(ω, source, material; kws...)

    # calculate eigenvectors
    MM_svd = svd(MM(k_eff))
    inds = findall(MM_svd.S .< tol)

    if isempty(inds)
        @warn("No eigenvectors found with the tolerance tol = $tol. Will use only one eigenvector with the eigenvalue $(MM_svd.S[end]), which should be less than tol.")
        inds = [length(MM_svd.S)]
    end

    #NOTE: MM(k_eff) ≈ MM_svd.U * diagm(0 => MM_svd.S) * MM_svd.Vt
    eigvectors = MM_svd.V[:,inds]

    # Reshape to separate different species and eigenvectors
    S = length(material.species)

    return reshape(eigvectors,(:,S,size(eigvectors,2)))
end
