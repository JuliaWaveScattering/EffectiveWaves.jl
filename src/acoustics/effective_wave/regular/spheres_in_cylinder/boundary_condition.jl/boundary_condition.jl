function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, source::AbstractSource{Acoustic{T,3}}, material::Material{Circle{T,2}}, ::TranslationSymmetry{3,T};
        basis_order::Int = 2,
        basis_field_order::Int = 2*basis_order,
        # source_basis_field_order::Int = basis_field_order,
        source_basis_field_order::Int = size(eigvectors)[3],
        kws...
    ) where T
    # source_basis_field_order is often chosen so that there is the same number of source coefficients a_n as the number of unknowns α_n
    # Before was: source_basis_field_order = min(basis_field_order,Int(round(sqrt(size(eigvectors)[end])))) - 1

    ρ = source.medium.ρ
    c = source.medium.c
    ρ0 = material.microstructure.medium.ρ
    c0 = material.microstructure.medium.c

    k = ω / c
    k0 = ω / c0

    species = material.microstructure.species
    S = length(species)
    rs = outer_radius.(species)

    # dim 1 is the (n,n1) indices, dim 2 is the species, dim 3 are the different eigenvectors
    # eigvectors = reshape(eigvectors,(:,S,size(eigvectors)[end]))

    R = outer_radius(material.shape)

    # Set the maximmum order the basis functions
    L = basis_order
    M = basis_field_order
    # L1 = Int(sqrt(size(eigvectors,1) / (L+1)^2) - 1)
    Minc = source_basis_field_order
    L2 = 2L + L1

    # the kernel use to wieight the species and the field's basis order.
    Ns = [
        (R - rs[j]) * kernelN3D(l1,k*(R - rs[j]), k_eff*(R - rs[j])) * number_density(species[j])
    for l1 = 0:L1, j in eachindex(species)] ./ (k^T(2) - k_eff^T(2))

    l1s = [l1 for l = 0:L for m = -l:l for l1 = 0:L1 for m1 = -l1:l1];

    vecs = [
        eigvectors[i] * Ns[l1s[i[1]]+1,i[2]]
    for i in CartesianIndices(eigvectors)];

    # sum over species
    vecs = sum(vecs, dims=2)
    vecs = reshape(vecs, size(eigvectors)[[1,3]])

    # Reflection and transmission coefficients
    γ = (ρ0 * c0) / (ρ * c)
    EW = EffectiveWaves

    Tran = [
        (γ*EW.diffhankelh1(m, k*R)*EW.besselj(m, k*R) - γ*EW.hankelh1(m, k*R)*EW.diffbesselj(m, k*R)) \
        (γ*EW.diffhankelh1(m, k*R)*EW.besselj(m, k0*R) - EW.hankelh1(m, k*R)*EW.diffbesselj(m, k0*R))
        for m = -Minc:Minc];

    Refl = [
        (EW.hankelh1(m, k*R)*EW.diffhankelh1(m, k0*R) - γ*EW.diffhankelh1(m, k*R))*EW.hankelh1(m, k0*R) \
        (γ*EW.diffhankelh1(m, k*R)*EW.besselj(m, k0*R) - γ*EW.hankelh1(m, k*R)*EW.diffbesselj(m, k0*R))
        for m = -Minc:Minc];

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
