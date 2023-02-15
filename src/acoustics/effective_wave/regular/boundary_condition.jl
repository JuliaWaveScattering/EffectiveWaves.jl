function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, source::AbstractSource{Acoustic{T,3}}, material::Material{Sphere{T,3}}, ::WithoutSymmetry{3};
        basis_order::Int = 2,
        basis_field_order::Int = 2*basis_order,
        # source_basis_field_order::Int = basis_field_order,
        source_basis_field_order::Int = Int(round(sqrt(size(eigvectors)[end]))) - 1,
        kws...
    ) where T
    # source_basis_field_order is often chosen so that there is the same number of source coefficients a_n as the number of unknowns α_n
    # Before was: source_basis_field_order = min(basis_field_order,Int(round(sqrt(size(eigvectors)[end])))) - 1

    if source.medium != material.microstructure.medium @error mismatched_medium end

    k = ω / source.medium.c

    species = material.microstructure.species
    S = length(species)
    rs = outer_radius.(species)

    # dim 1 is the (n,n1) indices, dim 2 is the species, dim 3 are the different eigenvectors
    # eigvectors = reshape(eigvectors,(:,S,size(eigvectors)[end]))

    R = outer_radius(material.shape)

    # Set the maximmum order the basis functions
    L = basis_order
    L1 = basis_field_order
    # L1 = Int(sqrt(size(eigvectors,1) / (L+1)^2) - 1)
    Linc = source_basis_field_order
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

function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, source::AbstractSource{Acoustic{T,3}}, material::Material{Sphere{T,3}}, ::AbstractAzimuthalSymmetry{3};
        basis_order::Int = 2,
        basis_field_order::Int = 4,
        # source_basis_field_order::Int = basis_field_order,
        source_basis_field_order::Int = size(eigvectors)[end] - 1,
        kws...
    ) where T
    # source_basis_field_order is often chosen so that there is the same number of source coefficients a_n as the number of unknowns α_n. Before was: min(basis_field_order, size(eigvectors)[end]) - 1

    k = ω / source.medium.c

    if source.medium != material.microstructure.medium @error mismatched_medium end
    species = material.microstructure.species
    S = length(species)
    rs = outer_radius.(species)

    # dim 1 is the (n,n1) indices, dim 2 is the species, dim 3 are the different eigenvectors
    # eigvectors = reshape(eigvectors,(:,S,size(eigvectors)[end]))

    R = outer_radius(material.shape)

    L = basis_order
    L1 = basis_field_order
    # L1 = Int((3 * size(eigvectors,1) - 3 - (2 - L) * L * (2 + L)) / (3 + 3*L*(2 + L)))
    Linc = source_basis_field_order

    # Linc = 2
    L2 = 2L + L1

    # the kernel use to wieight the species and the field's basis order.
    Ns = [
        (R - rs[j]) * kernelN3D(l1, k*(R - rs[j]), k_eff*(R - rs[j])) * number_density(species[j])
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

    # source = plane_source(psource.medium, psource.position, direction, psource.amplitude[1])
    # source.coefficients(Linc,zeros(3),ω)

    source_coefficients = regular_spherical_coefficients(source)(Linc,zeros(3),ω)

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

function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, source::AbstractSource{Acoustic{T,3}}, material::Material{Sphere{T,3}}, ::RadialSymmetry{3};
        basis_order::Int = 2,
        kws...
    ) where T

    k = ω / source.medium.c
    if source.medium != material.microstructure.medium @error mismatched_medium end

    species = material.microstructure.species
    S = length(species)
    rs = outer_radius.(species)

    R = outer_radius(material.shape)

    # the kernel use to weight the species and the field's basis order.
    F = sum(
        T(2 *(i[1] - 1) + 1) * (-one(T))^(i[1]-1) * eigvectors[i] * (R - rs[i[2]]) *
        kernelN3D(i[1] - 1, k*(R - rs[i[2]]), k_eff*(R - rs[i[2]])) * number_density(species[i[2]])
    for i in CartesianIndices(eigvectors)) / (k^T(2) - k_eff^T(2))

    # We expect there to only one component of a regular spherical wave expansion
    source_coefficients = regular_spherical_coefficients(source)(0,zeros(3),ω)

    forcing = source_coefficients / sqrt(T(4pi))

    α = F \ forcing

    return α
end

# Solve the boundary condition of the 2D radial symmetry case
function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, source::AbstractSource{Acoustic{T,2}}, material::Material{Sphere{T,2}}, ::RadialSymmetry{2};
        basis_order::Int = 2,
        kws...
    ) where T

    k = ω / source.medium.c
    species = material.microstructure.species
    S = length(species)
    rs = outer_radius.(species)

    R = outer_radius(material.shape)

    # the kernel use to weight the species and the field's basis order.
    F = 2pi*
    sum(
        kernelN2D(i[1]-1-basis_order, k*(R - rs[i[2]]), k_eff*(R - rs[i[2]])) *
        eigvectors[i] *
        number_density(species[i[2]])
    for i in CartesianIndices(eigvectors)) / (k^T(2) - k_eff^T(2))

    # We expect there to only one component of a regular cylindrical wave expansion
    # regular_spherical_coefficients(source)(n,ρ0,ω) computes g_n in g_nV_n(kρ-kρ0)
    source_coefficients = regular_spherical_coefficients(source)(0,zeros(2),ω)

    α = source_coefficients / F
    return α
end

# Solve the boundary condition of the spheres in cylinder case with 2 media
function solve_boundary_condition(ω::T, k_eff::Complex{T}, eigvectors::Array{Complex{T}}, source::AbstractSource{Acoustic{T,2}}, material::Material{Sphere{T,2}}, ::TranslationSymmetry{3,T};
        basis_order::Int = 3,
        basis_field_order::Int = 6,
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
        eigvectors[(lm_to_n(i[1],i[2]) - 1) * (2M + 1) + i[3] + M + 1, j, p] *
        Ns[i[3] + M1 + 1]
    for i in dn_m, j in eachindex(species), p = 1:(2Minc + 1)];

    # sum over species
    vecs = sum(vecs, dims=2)
    vecs = reshape(vecs, size(eigvectors)[[1,3]])

    # Reflection and transmission coefficients
    γ = (ρ0 * c0) / (ρ * c);

    Tran = [
        (γ*diffhankelh1(s, k*R)*besselj(s, k*R) - γ*hankelh1(s, k*R)*diffbesselj(s, k*R)) \
        (γ*diffhankelh1(s, k*R)*besselj(s, k0*R) - hankelh1(s, k*R)*diffbesselj(s, k0*R))
        for s = -Minc:Minc];

    Refl = [
        (hankelh1(s, k*R)*diffhankelh1(s, k0*R) - γ*diffhankelh1(s, k*R))*hankelh1(s, k0*R) \
        (γ*diffhankelh1(s, k*R)*besselj(s, k0*R) - hankelh1(s, k*R)*diffbesselj(s, k0*R))
        for s = -Minc:Minc];

    # Precomputation of spherical harmonic functions
    Ys = spherical_harmonics(L + M1, pi/2, 0.0)

    # Contributions from wall multiplescattering
    wall_reflections = [Refl[s + Minc + 1] *
        sum(l ->
            sum(m ->
                if abs(m - s) < l
                    Complex{T}(1im)^(l + m) * Ys[lm_to_n(l, m)] .*
                    vecs[(lm_to_n(l,m) - 1)*(2M1 + 1) + (m - s) + M1 + 1, p]
                else
                    Complex{T}(0.0)
                end
            ,-l:l)
        , 0:L) for s = -Minc:Minc, p = 1:(2Minc + 1)];

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
                        vecs[(lm_to_n(dl,dm) - 1)*(2M1 + 1) + (i[4] - m2) + M1 + 1, p]
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
        , 0:L) for i in n_n1];

    # Full matrix to be inverted
    matrix = (2pi^2/k0) * (particle_contribution + wall_contribution)

    # Reduction of matrix and forcing term
    square_matrix = zeros(Complex{T}, 2Minc+1, 2Minc+1)
    for q in 1:(2Minc + 1)
        for p in 1:(2Minc + 1)
            square_matrix[q, p] = sum(l ->
                sum(m ->
                    sum(dl ->
                        sum(dm ->
                            matrix[(lm_to_n(l,m) - 1)*(L+1)^2 + lm_to_n(dl,dm), p] * 
                            sum(ddl ->
                                gaunt_coefficient(ddl,q - Minc - 1,l,m,dl,dm)
                            , abs(q-Minc-1):Minc)
                        , -dl:dl)
                    , 0:L)
                , -l:l)
            , 0:L)
        end
    end

    reduced_forcing = zeros(Complex{T}, 2Minc + 1)
    for q in 1:(2Minc + 1)
        reduced_forcing[q] = sum(l ->
            sum(m ->
                sum(dl ->
                    sum(dm ->
                        forcing[(lm_to_n(l,m) - 1)*(L+1)^2 + lm_to_n(dl,dm)] *
                        sum(ddl ->
                            gaunt_coefficient(ddl,q - Minc - 1,l,m,dl,dm)
                        , abs(q - Minc - 1):Minc)
                    , -dl:dl)
                , 0:L)
            , -l:l)
        , 0:L)
    end

    # Computing normalization factors
    α = square_matrix \ reduced_forcing

    if rank(square_matrix) < 2Minc + 1
        @warn "Degeneracy encountered. The normalization factors computed in boundary conditions may not be reliable."
    end

    # Checking error in the solution of the system
    err = norm(reduced_forcing - square_matrix * α) / norm(reduced_forcing)

    if err > sqrt(eps(T))
        @warn "Extinction equation (like a boundary condition) was solved with an error: $err for the effective wavenumber: $k_eff"
    end

    return α

end