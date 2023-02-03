function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, source::AbstractSource{Acoustic{T,2}}, material::Material{Circle{T,2}}, ::TranslationSymmetry{3,T};
        basis_order::Int = 2,
        basis_field_order = 2 * basis_order,
        source_basis_field_order::Int = Int((size(eigvectors)[3] - 1) / 2),
        kws...
    ) where T
    # source_basis_field_order is often chosen so that there is the same number of source coefficients G_m as the number of unknowns α_n

    # Setting parameters
    ρ = source.medium.ρ
    ρ0 = material.microstructure.medium.ρ
    c = source.medium.c
    c0 = material.microstructure.medium.c
    k = ω / c
    k0 = ω / c0

    species = material.microstructure.species
    S = length(species)
    rs = outer_radius.(species)

    # dim 1 is the (n,n1) indices, dim 2 is the species, dim 3 are the different eigenvectors

    R = outer_radius(material.shape)

    # Set the maximmum order the basis functions
    L = basis_order
    M = basis_field_order
    Minc = source_basis_field_order
    M1 = L + Minc
    if M1 > M
        M1 = M
        @warn "Adjusting field order expansions for the system to be solvable"
    end

    # the kernel used to weight the species and the field's basis order.
    Ns = [
        kernelN2D(m, k0*(R - rs[j]), k_eff*(R - rs[j])) * number_density(species[j])
    for m = -M1:M1, j in eachindex(species)] ./ (k_eff^T(2) - k0^T(2))

    dn_m = [[dl,dm,m] for dl = 0:L for dm = -dl:dl for m = -M1:M1];
    lm_to_n = lm_to_spherical_harmonic_index

    vecs = [
        eigvectors[(lm_to_n(i[1],i[2]) - 1) * (2M1 + 1) + i[3] + M1 + 1, j, p] *
        Ns[i[3] + M1 + 1]
    for i in dn_m, j in eachindex(species), p = 1:(2Minc + 1)];

    # sum over species
    vecs = sum(vecs, dims=2)
    vecs = reshape(vecs, size(eigvectors)[[1,3]])

    # Reflection and transmission coefficients
    γ = (ρ0 * c0) / (ρ * c);

    Tran = [
        (γ*diffhankelh1(m, k*R)*besselj(m, k*R) - γ*hankelh1(m, k*R)*diffbesselj(m, k*R)) \
        (γ*diffhankelh1(m, k*R)*besselj(m, k0*R) - hankelh1(m, k*R)*diffbesselj(m, k0*R))
        for m = -Minc:Minc];

    Refl = [
        (hankelh1(m, k*R)*diffhankelh1(m, k0*R) - γ*diffhankelh1(m, k*R))*hankelh1(m, k0*R) \
        (γ*diffhankelh1(m, k*R)*besselj(m, k0*R) - hankelh1(m, k*R)*diffbesselj(m, k0*R))
        for m = -Minc:Minc];

    # Precomputation of spherical harmonic functions
    Ys = spherical_harmonics(L + M1, pi/2, 0.0)

    # Contributions from wall multiplescattering
    wall_reflections = [Refl[s + Minc + 1] *
        sum(l ->
            sum(m ->
                Complex{T}(1im)^(l + m) * Ys[lm_to_n(l, m)] .*
                vecs[(lm_to_n(l,m)-1)*(2M1+1) + (m-s) + M1 + 1, p]
            ,-l:l)
        , 0:(M1-Minc)) for s = -Minc:Minc, p = 1:(2Minc + 1)];

    n_n1 = [[l,m,l1,m1] for l = 0:L for m = -l:l for l1 = 0:L for m1 = -l1:l1];

    wall_contribution = [
        sum(dl ->
            sum(dm ->
                Complex{T}(1im)^(dl + dm) * Ys[lm_to_n(dl,dm)] *
                gaunt_coefficient(dl,dm,i[1],i[2],i[3],i[4]) *
                wall_reflections[dm + Minc + 1,p]
            , -dl:dl)
        , 0:L) for i in n_n1, p in 1:(2Minc + 1)];

    # Contributions from direct particle-particle multiplescattering
    particle_contribution = [
        Complex{T}(1im)^(-i[3] - i[4]) * Ys[lm_to_n(i[3],i[4])] *
        sum(l2 ->
            sum(m2 ->
                Complex{T}(1im)^(l2 + m2) * Ys[lm_to_n(l2,m2)] *
                sum(dl ->
                    sum(dm ->
                        gaunt_coefficient(dl,dm,i[1],i[2],l2,m2) *
                        vecs[(lm_to_n(dl,dm)-1)*(2M+1) + (i[4]-m2) + M1 + 1, p]
                    , -dl:dl)
                , 0:L)
            , -l2:l2)
        , 0:L) for i in n_n1, p in 1:(2Minc + 1)];

    # Incident wave coefficients
    source_coefficients = Tran .* regular_spherical_coefficients(source)(Minc,zeros(2),ω)

    # Forcing term
    forcing = [Complex{T}(-1) *
        sum(dl ->
            sum(dm ->
                Complex{T}(1im)^(dl + dm) * Ys[lm_to_n(dl,dm)] *
                gaunt_coefficient(dl,dm,i[1],i[2],i[3],i[4]) *
                source_coefficients[dm + Minc + 1]
            , -dl:dl)
        , 0:L)
        for i in n_n1]

    # Full matrix to be inverted
    matrix = (2pi^2/k0) * (particle_contribution + wall_contribution)

    # Reduction of matrix and forcing term
    square_matrix = [
        sum(l ->
            sum(m ->
                sum(dl ->
                    sum(dm ->
                        matrix[(lm_to_n(l,m)-1)*(L+1)^2 + lm_to_n(dl,dm), p] *
                        sum(ddl ->
                            gaunt_coefficient(ddl,ddm,l,m,dl,dm)
                        , abs(ddm):Minc)
                    , -dl:dl)
                , 0:L)
            , -l:l)
        , 0:L) for ddm in -Minc:Minc, p in 1:(2Minc + 1)]

    reduced_forcing = [
        sum(l ->
            sum(m ->
                sum(dl ->
                    sum(dm ->
                        forcing[(lm_to_n(l,m)-1)*(L+1)^2 + lm_to_n(dl,dm)] *
                        sum(ddl ->
                            gaunt_coefficient(ddl,ddm,l,m,dl,dm)
                        , abs(ddm):Minc)
                    , -dl:dl)
                , 0:L)
            , -l:l)
        , 0:L) for ddm in -Minc:Minc]

    # Computing normalization factors
    α = square_matrix \ reduced_forcing

    # Checking error in the solution of the system
    err = norm(reduced_forcing - square_matrix * α) / norm(reduced_forcing)

    if err > sqrt(eps(T))
        @warn "Extinction equation (like a boundary condition) was solved with an error: $err for the effective wavenumber: $k_eff"
    end

    return α

end
