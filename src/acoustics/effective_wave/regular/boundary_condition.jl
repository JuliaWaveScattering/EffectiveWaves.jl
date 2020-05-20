function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, source::Source{T,Acoustic{T,3}}, material::Material{3,Sphere{T}};
        basis_order::Int = 2,
        kws...
    ) where T

    k = real(ω / source.medium.c)

    species = material.species
    S = length(species)
    as = outer_radius.(species)

    # dim 1 is the (n,n1) indices, dim 2 is the species, dim 3 are the different eigenvectors
    eigvectors = reshape(eigvectors,(:,S,size(eigvectors,2)))

    R = outer_radius(material.shape)

    L = basis_order
    L1 = Int(sqrt(size(eigvectors,1) / (L+1)^2) - 1)

    # setting Linc = P guarantees that there could be a one-to-one relationship between the incident wave and the solution.
    Linc = Int(sqrt(size(eigvectors)[end]) - 1)
    L2 = 2L + L1

    # the kernel use to wieight the species and the field's basis order.
    Ns = [
        (R - as[j]) * kernelN3D(l1,k*(R - as[j]), k_eff*(R - as[j])) * number_density(species[j])
    for l1 = 0:L1, j in eachindex(species)] ./ (k^T(2) - k_eff^T(2))

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

    source_coefficients = source.coefficients(Linc,zeros(3),ω)

    forcing = [
        - sum(
            [gaunt_coefficient(dl,dm,l,m,l2,m2) for dl = 0:Linc for dm = -dl:dl] .*
            source_coefficients
        )
    for l = 0:L for m = -l:l
    for l2 = 0:L2 for m2 = -l2:l2]

    α = (extinction_matrix * eigvectors) \ forcing

    err = norm(forcing - extinction_matrix * eigvectors * α) / norm(forcing)

    if err > sqrt(eps(T))
        error("Extinction equation (like a boundary condition) was solved with an error: $err ")
    end

    return α

end

function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, psource::PlaneSource{T,3,1,Acoustic{T,3}}, material::Material{3,Sphere{T}};
        basis_order::Int = 2,
        kws...
    ) where T

    k = real(ω / psource.medium.c)

    species = material.species
    S = length(species)
    as = outer_radius.(species)

    # dim 1 is the (n,n1) indices, dim 2 is the species, dim 3 are the different eigenvectors
    eigvectors = reshape(eigvectors,(:,S,size(eigvectors,2)))

    R = outer_radius(material.shape)

    L = basis_order
    L1 = Int((3 * size(eigvectors,1) - 3 - (2 - L) * L * (2 + L)) / (3 + 3*L*(2 + L)))

    # Below guarantees that there is the same number of possible a_n as α_n
    Linc = size(eigvectors,3) - 1
    L2 = 2L + L1

    # the kernel use to wieight the species and the field's basis order.
    Ns = [
        (R - as[j]) * kernelN3D(l1,k*(R - as[j]), k_eff*(R - as[j])) * number_density(species[j])
    for l1 = 0:L1, j in eachindex(species)] ./ (k^T(2) - k_eff^T(2))

    l1s = [l1 for l = 0:L for m = -l:l for l1 = abs(m):L1];

    eigvectors = [
        eigvectors[i] * Ns[l1s[i[1]]+1,i[2]]
    for i in CartesianIndices(eigvectors)];

    # sum over species
    eigvectors = sum(eigvectors, dims=2)
    eigvectors = reshape(eigvectors, size(eigvectors)[[1,3]])

    # The extinction_matrix
    function gaunt2(dl,dm,l1,l,m,l2)::Complex{T}
        minl3 = max(abs(dl-l),abs(l2-l1),abs(dm-m))
        maxl3 = min(abs(dl+l),abs(l2+l1))
        return if minl3 <= maxl3
            - sum(l3 ->
                gaunt_coefficient(dl,dm,l,m,l3,dm-m) *
                gaunt_coefficient(l2,-m,l1,-dm,l3,dm-m)
            , minl3:maxl3)
        else
            zero(Complex{T})
        end
    end

    # in this form: extinction_matrix[(n,n2),(dn,n1)] ==  gaunt2(dl,dm,l1,m1,l,m,l2,m2)
    extinction_matrix = [
            gaunt2(dl,dm,l1,l,m,l2)
        for dl = 0:L for dm = -dl:dl
        for l1 = abs(dm):L1
    for l = 0:L for m = -l:l
    for l2 = abs(m):L2]

    len = Int(1 - L*(2 + L)*(L - 3*L1 - 2)/3 + L1)
    extinction_matrix = reshape(extinction_matrix, (:,len))

    # we expect the source direction to be aligned with the z-axis
    direction = SVector(zero(T),zero(T),one(T))

    if norm(abs.(psource.direction) - abs.(direction)) > eps(T)
        warn("Plane wave source is not aligned with the z-axis. We return results with an axis that makes the source direction aligned with the z-axis")
    end

    source = plane_source(psource.medium, psource.position, direction, psource.amplitude[1])

    source_coefficients = source.coefficients(Linc,zeros(3),ω)

    forcing = [
        - sum(
            [gaunt_coefficient(dl,dm,l,m,l2,-m) for dl = 0:Linc for dm = -dl:dl] .*
            source_coefficients
        )
    for l = 0:L for m = -l:l
    for l2 = abs(m):L2]

    # If the forcing was exactly a plane-source in the z-axis direction:
    # forcing2 = - [
    #     (m != 0 || m2 != 0) ? 0.0im : 4pi * (T(1)*im)^(l2+l) * sqrt((2l2+1)*(2l+1))
    # for l = 0:L for m = -l:l
    # for l2 = 0:L2 for m2 = -l2:l2]
    #
    # findall( abs.(forcing) .> 1e-4)
    # norm(forcing - forcing2) / norm(forcing2)

    α = (extinction_matrix * eigvectors) \ forcing

    err = norm(forcing - extinction_matrix * eigvectors * α) / norm(forcing)

    if err > sqrt(eps(T))
        error("Extinction equation (like a boundary condition) was solved with an error: $err ")
    end

    return α

end
