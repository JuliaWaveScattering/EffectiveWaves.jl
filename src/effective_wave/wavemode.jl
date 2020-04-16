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

    extinction_matrix, forcing = boundary_condition_system(ω, k_eff, source, material; kws...)
    # extinction_matrix, forcing = boundary_condition_system(ω, k_eff, source, material; basis_order = basis_order)

    α = (extinction_matrix * eigvectors) \ forcing

    # normalise the eigvectors such that their sum satisfies the boundary condition
    eigvectors = [eigvectors[i] * α[i[2]] for i in CartesianIndices(eigvectors)]

    # Reshape to separate different species and eigenvectors
    S = length(material.species)

    return reshape(eigvectors,(:,S,size(eigvectors,2)))
end


# "The average effective transmitted scattering amplitudes for a given effective wavenumber k_eff. Assumes there exists only one k_eff.
# The function returns an array A, where
# AA(x,y,m,s) = im^m*exp(-im*m*θ_eff)*A[m + max_basis_order +1,s]*exp(im*k_eff*(cos(θ_eff)*x + sin(θin)*y))
# where (x,y) are coordinates in the halfspace, m-th hankel order, s-th species,  and AA is the ensemble average scattering coefficient."
# function mode_amplitudes(ω::T, k_eff::Complex{T}, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, material::Material{Dim,Halfspace{T,Dim}};
#         tol::T = 1e-5,
#         kws...) where {T<:Number,Dim}
#
#     MM = eigensystem(ω, psource, material; kws...)
#
#     # calculate effective amplitudes
#     MM_svd = svd(MM(k_eff))
#     if MM_svd.S[end] > tol*T(100)
#         @warn("The effective wavenumber gave an eigenvalue $(MM_svd.S[end]), which would be zero if the solution was exact.")
#     end
#
#     S = length(material.species)
#
#     A_null = MM_svd.V[:,end] # eignvector of smallest eigenvalue
#     A_null = reshape(A_null, (:,S)) # A_null[:,j] = [A[-ho,j],A[-ho+1,j],...,A[ho,j]]
#
#     return A_null
# end
