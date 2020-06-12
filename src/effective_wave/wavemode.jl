eigensystem(ω::T, source::AbstractSource{T}, material::Material; kws...) where T<:AbstractFloat = eigensystem(ω, source.medium, material.species, setupsymmetry(source,material); kws...)

"Calculates the effective wavenumbers and return Vector{EffectivePlaneWaveMode}."
function WaveModes(ω::T, source::AbstractSource, material::Material{Dim,S,Sps}; kws...) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}} # without the parametric types we get a "Unreachable reached" error

    k_effs = wavenumbers(ω, source.medium, material.species; kws... )
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

    vecs = eigenvectors(ω, wavenumber, source, material; kws...)

    return EffectiveRegularWaveMode(ω, wavenumber, source, material, vecs; kws...)
end

function WaveMode(ω::T, wavenumber::Complex{T}, psource::PlaneSource{T,Dim,1}, material::Material{Dim,Halfspace{T,Dim}};
    tol::T = 1e-6, kws...) where {T,Dim}

    direction = transmission_direction(wavenumber, (ω / psource.medium.c) * psource.direction, material.shape.normal; tol = tol)

    vecs = eigenvectors(ω, wavenumber, psource, material; kws...)

    return EffectivePlaneWaveMode(ω, wavenumber, direction, vecs)
end

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

    α = solve_boundary_condition(ω, k_eff, eigvectors, source, material;
            kws...
    )

    # After this normalisation, sum(eigvectors, dims=2) will satisfy the boundary conditions
    eigvectors = [eigvectors[i] * α[i[2]] for i in CartesianIndices(eigvectors)]

    # Reshape to separate different species and eigenvectors
    S = length(material.species)

    return reshape(eigvectors,(:,S,size(eigvectors,2)))
end
