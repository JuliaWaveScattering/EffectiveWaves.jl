function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, source::Source{T,Acoustic{T,3}}, material::Material{3,Sphere{T,3}};
        basis_order::Int = 2,
        basis_field_order::Int = 4,
        source_basis_field_order::Int = Int(round(sqrt(size(eigvectors)[end]))) - 1,
        kws...
    ) where T
    # source_basis_field_order is often chosen so that there is the same number of source coefficients a_n as the number of unknowns α_n
    # Before was: source_basis_field_order = min(basis_field_order,Int(round(sqrt(size(eigvectors)[end])))) - 1

    k = real(ω / source.medium.c)

    species = material.species
    S = length(species)
    as = outer_radius.(species)

    # dim 1 is the (n,n1) indices, dim 2 is the species, dim 3 are the different eigenvectors
    eigvectors = reshape(eigvectors,(:,S,size(eigvectors)[end]))

    R = outer_radius(material.shape)

    # Set the maximmum order the basis functions
    L = basis_order
    L1 = basis_field_order
    # L1 = Int(sqrt(size(eigvectors,1) / (L+1)^2) - 1)
    Linc = source_basis_field_order
    L2 = 2L + L1

    # the kernel use to wieight the species and the field's basis order.
    Ns = [
        (R - as[j]) * kernelN3D(l1,k*(R - as[j]), k_eff*(R - as[j])) * number_density(species[j])
    for l1 = 0:L1, j in eachindex(species)] ./ (k^T(2) - k_eff^T(2))

    l1s = [l1 for l = 0:L for m = -l:l for l1 = 0:L1 for m1 = -l1:l1];

    vecs = [
        eigvectors[i] * Ns[l1s[i[1]]+1,i[2]]
    for i in CartesianIndices(eigvectors)];

    # sum over species
    vecs = sum(vecs, dims=2)
    vecs = reshape(vecs, size(eigvectors)[[1,3]])

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

    source_coefficients = regular_spherical_coefficients(source)(Linc,zeros(3),ω)

    forcing = [
        - sum(
            [gaunt_coefficient(dl,dm,l,m,l2,m2) for dl = 0:Linc for dm = -dl:dl] .*
            source_coefficients
        )
    for l = 0:L for m = -l:l
    for l2 = 0:L2 for m2 = -l2:l2]

    α = (extinction_matrix * vecs) \ forcing

    err = norm(forcing - extinction_matrix * vecs * α) / norm(forcing)

    if err > sqrt(eps(T))
        @warn "Extinction equation (like a boundary condition) was solved with an error: $err for the effective wavenumber: $k_eff"
    end

    return α

end

function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, psource::PlaneSource{T,3,1,Acoustic{T,3}}, material::Material{3,Sphere{T,3}};
        basis_order::Int = 2,
        basis_field_order::Int = 4,
        source_basis_field_order::Int = size(eigvectors)[end] - 1,
        kws...
    ) where T
    # source_basis_field_order is often chosen so that there is the same number of source coefficients a_n as the number of unknowns α_n. Before was: min(basis_field_order, size(eigvectors)[end]) - 1

    k = real(ω / psource.medium.c)

    species = material.species
    S = length(species)
    as = outer_radius.(species)

    # dim 1 is the (n,n1) indices, dim 2 is the species, dim 3 are the different eigenvectors
    eigvectors = reshape(eigvectors,(:,S,size(eigvectors,2)))

    R = outer_radius(material.shape)

    L = basis_order
    L1 = basis_field_order
    # L1 = Int((3 * size(eigvectors,1) - 3 - (2 - L) * L * (2 + L)) / (3 + 3*L*(2 + L)))
    Linc = source_basis_field_order

    # Linc = 2
    L2 = 2L + L1

    # the kernel use to wieight the species and the field's basis order.
    Ns = [
        (R - as[j]) * kernelN3D(l1,k*(R - as[j]), k_eff*(R - as[j])) * number_density(species[j])
    for l1 = 0:L1, j in eachindex(species)] ./ (k^T(2) - k_eff^T(2))

    l1s = [l1 for l = 0:L for m = -l:l for l1 = abs(m):L1];

    vecs = [
        eigvectors[i] * Ns[l1s[i[1]]+1,i[2]]
    for i in CartesianIndices(eigvectors)];

    # sum over species
    vecs = sum(vecs, dims=2)
    vecs = reshape(vecs, size(vecs)[[1,3]])

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

    extinction_matrix = reshape(extinction_matrix, (:,size(vecs,1)))

    # we expect the source direction to be aligned with the z-axis
    direction = SVector(zero(T),zero(T),one(T))

    if norm(abs.(psource.direction) - abs.(direction)) > eps(T)
        warn("Plane wave source is not aligned with the z-axis. Will return results with an axis that makes the source direction aligned with the z-axis")
    end

    # source = plane_source(psource.medium, psource.position, direction, psource.amplitude[1])
    # source.coefficients(Linc,zeros(3),ω)

    source_coefficients = regular_spherical_coefficients(psource)(Linc,zeros(3),ω)

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

    α = (extinction_matrix * vecs) \ forcing

    err = norm(forcing - extinction_matrix * vecs * α) / norm(forcing)

    if err > sqrt(eps(T))
        @warn "Extinction equation (like a boundary condition) was solved with an error: $err for the effective wavenumber: $k_eff"
    end

    return α

end
