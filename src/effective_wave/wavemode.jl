eigensystem(ω::T, source::AbstractSource{T}, material::Material; kws...) where T<:AbstractFloat = eigensystem(ω, source.medium, material.species, setupsymmetry(source,material); kws...)

"Calculates the effective wavenumbers and return Vector{EffectivePlaneWaveMode}."
function wavemodes(ω::T, source::AbstractSource, material::Material{Dim,S,Sps}; kws...) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}} # without the parametric types we get a "Unreachable reached" error

    k_effs = wavenumbers(ω, source.medium, material.species; kws... )
    wave_effs = [
        wavemode(ω, k_eff, source, material; kws...)
    for k_eff in k_effs]

    return wave_effs
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

    α = solve_boundary_condition(ω, k_eff, eigvectors, source, material; kws...)

    # The sum of these normalised vectors will now satisfy the boundary conditions
    eigvectors = [eigvectors[i] * α[i[2]] for i in CartesianIndices(eigvectors)]

    # Reshape to separate different species and eigenvectors
    S = length(material.species)

    return reshape(eigvectors,(:,S,size(eigvectors,2)))
end
