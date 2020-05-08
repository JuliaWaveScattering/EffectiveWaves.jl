function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, source::Source{T,Acoustic{T,3}}, material::Material{3,Sphere{T}};
        basis_order::Int = 2,
        basis_field_order::Int = 4,
        kws...
    ) where {T<:Number,Dim}

    k = real(ω / source.medium.c)

    species = material.species
    S = length(species)
    as = outer_radius.(species)

    R = outer_radius(material.shape)

    L = basis_order
    L1 = basis_field_order

    Linc = estimate_regular_basisorder(typeof(source.medium), R * k )
    L2 = L + Linc
    # L2 = L1
    # or is it this:

    # the kernel use to wieight the species and the field's basis order.
    Ns = [
        (R - as[j]) * kernelN3D(l1,k*(R - as[j]), k_eff*(R - as[j])) * number_density(species[j])
    for l1 = 0:basis_field_order, j in eachindex(species)] ./ (k^T(2) - k_eff^T(2))

    # dim 1 is the (n,n1) indices, dim 2 is the species, dim 3 are the different eigenvectors
    eigvectors = reshape(eigvectors,(:,S,size(eigvectors,2)))

    l1s = [l1 for l = 0:L for m = -l:l for l1 = 0:L1 for m1 = -l1:l1];

    eigvectors = [
        eigvectors[i] * Ns[l1s[i[1]]+1,i[2]]
    for i in CartesianIndices(eigvectors)];

    # sum over species
    eigvectors = sum(eigvectors, dims=2)
    eigvectors = reshape(eigvectors, size(eigvectors)[[1,3]])

    # The extinction_matrix
    function gaunt2(dl,dm,l1,m1,l,m,l2,m2)::Complex{T}
        minl3 = max(abs(dl-l),abs(l2-l1),abs(dm-m))
        maxl3 = min(abs(dl+l),abs(l2+l1))
        return if minl3 <= maxl3
            - sum(l3 ->
                gaunt_coefficient(dl,dm,l,m,l3,dm-m) *
                gaunt_coefficient(l2,m2,l1,m1,l3,dm-m)
            , minl3:maxl3)
        else
            zero(Complex{T})
        end
    end

    # in this form: extinction_matrix[(n,n2),(dn,n1)] ==  gaunt2(dl,dm,l1,m1,l,m,l2,m2)
    extinction_matrix = [
            gaunt2(dl,dm,l1,m1,l,m,l2,m2)
        for dl = 0:L for dm = -dl:dl
        for l1 = 0:L1 for m1 = -l1:l1
    for l = 0:L for m = -l:l
    for l2 = 0:L2 for m2 = -l2:l2]

    len = (L1+1)^2 * (L+1)^2
    extinction_matrix = reshape(extinction_matrix, (:,len))

    # Linc = max(L2+2,Linc)

    source_coefficients = source.coefficients(Linc,zeros(3),ω)

    forcing = [
        - sum(
            [gaunt_coefficient(dl,dm,l,m,l2,m2) for dl = 0:Linc for dm = -dl:dl] .*
            source_coefficients
        )
    for l = 0:L for m = -l:l
    for l2 = 0:L2 for m2 = -l2:l2]

    forcing2 = - [
        (m != 0 || m2 != 0) ? 0.0im : 4pi * (T(1)*im)^(l2+l) * sqrt((2l2+1)*(2l+1))
    for l = 0:L for m = -l:l
    for l2 = 0:L2 for m2 = -l2:l2]

    norm(forcing - forcing2) / norm(forcing2)

    forcing = forcing2

    α = (extinction_matrix * eigvectors) \ forcing

    norm(forcing - extinction_matrix * eigvectors * α) / norm(forcing)
    inds = findall( abs.(forcing - extinction_matrix * eigvectors * α) .> 1e-4)
    forcing[inds]

    length(inds)
    # 17 for L1 = 8
    # 31, 29 for L1 = 7
    # 27 for L1 = 4
    # 14 for L1 = 2


    return α

end

# function wavemode(ω::T, k_eff::Complex{T}, source::Source{T,Acoustic{T,Dim}}, material::Material{Dim,Sphere{T}};
#         tol::T = 1e-6, kws...
#     ) where {T<:AbstractFloat,Dim}
#
#     k = ω/psource.medium.c
#
#     direction = transmission_direction(k_eff, (ω / psource.medium.c) * psource.direction, material.shape.normal; tol = tol)
#
#     amps = eigenvectors(ω, k_eff, psource, material; tol= tol)
#
#     return EffectiveRegularWaveMode(ω, k_eff, direction, amps)
# end
