kernelM(l::Int,x, y) = x * diffsbesselj(l,x) * sbesselj(l,y) - y * sbesselj(l,x) * diffsbesselj(l,y)

function scattering_field(wavemode::EffectiveRegularWaveMode{T,Acoustic{T,3},WithoutSymmetry{3}}) where T
    L = wavemode.basis_order
    L1 = wavemode.basis_field_order

    i2s = 1:size(wavemode.eigenvectors,2)
    i3s = 1:size(wavemode.eigenvectors,3)

    nn1_to_n = [
            basisorder_to_basislength(Acoustic{T,3},l-1) + m + l + 1
    for l = 0:L for m = -l:l for l1 = 0:L1 for m1 = -l1:l1];

    nn1_to_n1 = [
            basisorder_to_basislength(Acoustic{T,3},l1-1) + m1 + l1 + 1
    for l = 0:L for m = -l:l for l1 = 0:L1 for m1 = -l1:l1];

    indices_same_n = collect(groupby(nn1 -> nn1_to_n[nn1] , 1:length(nn1_to_n)))

    v = regular_basis_function(wavemode.wavenumber, wavemode.source.medium)

    function scat_field(x::AbstractVector{T})
        vs =  v(L1, x)

        vecs = [
            sum(wavemode.eigenvectors[ind,i2,i3] .* vs[nn1_to_n1[ind]]) # note here:  vs[nn1_to_n1[ind]]] == vs
        for ind in indices_same_n, i2 in i2s, i3 in i3s];

        # Sum over different eigvectors
        vecs = sum(vecs, dims = 3)[:,:];

        # for only 1 species the code below gives the same
        # nn1_to_n1 = [
        #         basisorder_to_basislength(Acoustic{T,3},l1-1) + m1 + l1 + 1
        # for l = 0:L for m = -l:l for l1 = 0:L1 for m1 = -l1:l1];
        #
        # vecs2 = [
        #     wavemode.eigenvectors[i[1],i[2],i[3]] * vs[nn1_to_n1[i[1]]]
        # for i in CartesianIndices(wavemode.eigenvectors)];
        #
        # vecs2 = sum(vecs2, dims = 3);
        # vecs2 = [sum(vecs2[inds]) for inds in indices_same_n]
        # vecs ≈ vecs2
        return vecs
    end

end

function scattering_field(wavemode::EffectiveRegularWaveMode{T,Acoustic{T,3},Sym}) where {T, Sym<:AbstractAzimuthalSymmetry{3}}
    L = wavemode.basis_order
    L1 = wavemode.basis_field_order

    i2s = 1:size(wavemode.eigenvectors,2)
    i3s = 1:size(wavemode.eigenvectors,3)

    nl1_to_n = [
            basisorder_to_basislength(Acoustic{T,3},l-1) + m + l + 1
    for l = 0:L for m = -l:l for l1 = abs(m):L1];

    nl1_to_n1 = [
            basisorder_to_basislength(Acoustic{T,3},l1-1) - m + l1 + 1
    for l = 0:L for m = -l:l for l1 = abs(m):L1];

    indices_same_n = collect(groupby(nl1 -> nl1_to_n[nl1] , 1:length(nl1_to_n)))

    v = regular_basis_function(wavemode.wavenumber, wavemode.source.medium)

    function scat_field(x::AbstractVector{T})
        vs =  v(L1, x)

        vecs = [
            sum(wavemode.eigenvectors[ind,i2,i3] .* vs[nl1_to_n1[ind]])
        for ind in indices_same_n, i2 in i2s, i3 in i3s];

        # Sum over different eigvectors
        vecs = sum(vecs, dims = 3)[:,:];

        return vecs
    end

end

"""
    scattering_field(wavemode::EffectiveRegularWaveMode{T,Acoustic{T,3},RadialSymmetry{3}})

In general ``F_{nn'} = \\delta_{m',-m} \\delta_{\\ell,\\ell'}(-1)^m F_{(\\ell,0),(\\ell,0)}``, so the scattering field is given by

``\\langle f_n \\rangle (x) = \\sum_{n'} F_{nn'} \\mathrm v_{n'}(k_* x) =  F_{(\\ell,0),(\\ell,0)} (-1)^m \\mathrm v_{(\\ell,-m)}(k_* x)``
"""
function scattering_field(wavemode::EffectiveRegularWaveMode{T,Acoustic{T,3},RadialSymmetry{3}}) where T
    L = wavemode.basis_order
    len = basisorder_to_basislength(Acoustic{T,3},L)

    ls, ms = spherical_harmonics_indices(L)
    inds = lm_to_spherical_harmonic_index.(ls,-ms)

    v = regular_basis_function(wavemode.wavenumber, wavemode.source.medium)

    function scat_field(x::AbstractVector{T})
        vs =  v(L, x)
        vs = vs[inds] .* (-one(T)) .^ ms

        vecs = [
            sum(wavemode.eigenvectors[ls[n]+1,i2,:]) * vs[n]
        for n = 1:len, i2 = 1:size(wavemode.eigenvectors,2)];

        return vecs
    end
end



function material_scattering_coefficients(wavemode::EffectiveRegularWaveMode{T,Acoustic{T,3},WithoutSymmetry{3}}) where T

    # Unpacking parameters
    k = wavemode.ω / wavemode.source.medium.c

    species = wavemode.material.microstructure.species
    S = length(species)
    rs = outer_radius.(species)

    R = outer_radius(wavemode.material.shape)
    k_eff = wavemode.wavenumber

    L = wavemode.basis_order
    L1 = wavemode.basis_field_order
    dL = L + L1

    # Pre-calculating
    Ml1s = [
        (R - rs[j]) * kernelM(l1, k*(R - rs[j]), k_eff*(R - rs[j])) * number_density(species[j])
    for l1 = 0:L1, j in eachindex(species)] ./ (k_eff^T(2) - k^T(2))

    nn1_to_l1 = [l1 for l = 0:L for m = -l:l for l1 = 0:L1 for m1 = -l1:l1];

    # Sum over species and different eigvectors
    vecs = [
        wavemode.eigenvectors[i] * Ml1s[nn1_to_l1[i[1]]+1,i[2]]
    for i in CartesianIndices(wavemode.eigenvectors)]

    vec = reshape(sum(vecs, dims = (2,3)),:)

    # Sum over wave basis components
    gaunt = [
            gaunt_coefficient(dl,dm,l,m,l1,m1)
        for l = 0:L for m = -l:l
        for l1 = 0:L1 for m1 = -l1:l1
    for dl = 0:dL for dm = -dl:dl]

    # gaunt[(dn),(n,l1)] = gaunt_coefficient(dl,dm,l,m,l1,m1)
    gaunt = reshape(gaunt, (:,length(vec)))

    return gaunt * vec
end

function material_scattering_coefficients(wavemode::EffectiveRegularWaveMode{T,Acoustic{T,3},Sym}) where {T,Sym<:AbstractAzimuthalSymmetry{3}}

    # Unpacking parameters
    k = wavemode.ω / wavemode.source.medium.c

    species = wavemode.material.microstructure.species
    S = length(species)
    rs = outer_radius.(species)

    R = outer_radius(wavemode.material.shape)
    k_eff = wavemode.wavenumber

    L = wavemode.basis_order
    L1 = wavemode.basis_field_order
    dL = L + L1

    # Pre-calculating
    Ml1s = [
        (R - rs[j]) * kernelM(l1, k*(R - rs[j]), k_eff*(R - rs[j])) * number_density(species[j])
    for l1 = 0:L1, j in eachindex(species)] ./ (k_eff^T(2) - k^T(2))

    nl1_to_l1 = [l1 for l = 0:L for m = -l:l for l1 = abs(m):L1];

    # Sum over species and different eigvectors
    vecs = [
        wavemode.eigenvectors[i] * Ml1s[nl1_to_l1[i[1]]+1,i[2]]
    for i in CartesianIndices(wavemode.eigenvectors)]

    vec = reshape(sum(vecs, dims = (2,3)),:)

    # Sum over wave basis components
    gaunt = [
            gaunt_coefficient(dl,dm,l,m,l1,-m)
        for l = 0:L for m = -l:l
        for l1 = abs(m):L1
    for dl = 0:dL for dm = -dl:dl]

    gaunt = reshape(gaunt, (:,length(vec)))

    return gaunt * vec
end

function material_scattering_coefficients(wavemode::EffectiveRegularWaveMode{T,Acoustic{T,3},RadialSymmetry{3}}) where T

    # Unpacking parameters
    k = wavemode.ω / wavemode.source.medium.c
    k_eff = wavemode.wavenumber
    R = outer_radius(wavemode.material.shape)

    species = wavemode.material.microstructure.species
    S = length(species)
    rs = outer_radius.(species)

    L = wavemode.basis_order

    Fscat0 = sum(
        sqrt(T(4pi)) * wavemode.eigenvectors[i] * (2*(i[1]-1) + 1) * T(-1)^(i[1]-1) *
        (R - rs[i[2]]) * kernelM(i[1]-1, k*(R - rs[i[2]]), k_eff*(R - rs[i[2]])) * number_density(species[i[2]])
    for i in CartesianIndices(wavemode.eigenvectors)) / (k_eff^T(2) - k^T(2))

    return [Fscat0]
end

# Compute the coefficients Fs of the averaged scattered field in the case of spheres in a cylinder with two mwdia
function material_scattering_coefficients(wavemode::EffectiveRegularWaveMode{T,Acoustic{T,2},TranslationSymmetry{3,T}};
        only_particle_contribution::Bool = false,
        kws...
    ) where T

    # Unpacking parameters
    ρ = wavemode.source.medium.ρ
    ρ0 = wavemode.material.microstructure.medium.ρ
    c = wavemode.source.medium.c
    c0 = wavemode.material.microstructure.medium.c
    k = wavemode.ω / c
    k0 = wavemode.ω / c0
    k_eff = wavemode.wavenumber
    R = outer_radius(wavemode.material.shape)

    species = wavemode.material.microstructure.species
    Sp = length(species)
    rs = outer_radius.(species)

    eigenvectors = wavemode.eigenvectors
    P = size(eigenvectors)[3]

    L = wavemode.basis_order
    M1 = wavemode.basis_field_order
    Minc = M1 - L
    if Minc <= 0
        @warn "Not enough terms in field expansion, please get a higher value for basis_field_order = $basis_field_order."
    end

    # Precompliling spherical harmonics
    Ys = spherical_harmonics(L, pi/2, 0.0)
    lm_to_n = lm_to_spherical_harmonic_index

    # Computing transmission and reflection coefficients
    γ = (ρ0 * c0) / (ρ * c);

    Tran = [
        (diffhankelh1(s, k0*R)*besselj(s, k0*R) - hankelh1(s, k0*R)*diffbesselj(s, k0*R)) \
        (γ*diffhankelh1(s, k*R)*besselj(s, k0*R) - hankelh1(s, k*R)*diffbesselj(s, k0*R))
        for s = -Minc:Minc];

    Refl = [
        (besselj(s, k*R)*diffbesselj(s, k0*R) - γ*diffbesselj(s, k*R))*besselj(s, k0*R) \
        (γ*diffhankelh1(s, k*R)*besselj(s, k0*R) - hankelh1(s, k*R)*diffbesselj(s, k0*R))
        for s = -Minc:Minc];

    # Computing internal field contributions
    particle_contribution = zeros(Complex{T}, 2Minc + 1)
    for s in -Minc:Minc
        particle_contribution[s + Minc + 1] = Complex{T}(2*pi^2/(k0 * (k_eff^2 - k0^2))) *
        sum(l ->
            sum(m ->
                Complex{T}(1im)^(m - l) * Ys[lm_to_n(l, m)] *
                sum(i ->
                    number_density(species[i]) * kernelN2D(s - m, k0*(R - rs[i]), k_eff*(R - rs[i])) *
                    sum(p ->
                        eigenvectors[(lm_to_n(l, m) - 1)*(2M1 + 1) + (s - m) + M1 + 1, i, p]
                    , 1:P)
                , 1:Sp)
            , -l:l)
        , 0:L)
    end

    # Incident wave contributions
    if only_particle_contribution
        wall_contribution = zeros(Complex{T}, 2Minc + 1)
    else
        wall_contribution = regular_spherical_coefficients(wavemode.source)(Minc,zeros(2),wavemode.ω)
    end

    Fscat = (Tran .* particle_contribution) + (Refl .* wall_contribution)

    return Fscat
end

# Compute the coefficient F0 of the averaged scattered field in the case of the 2D radial symmetry case
function material_scattering_coefficients(wavemode::EffectiveRegularWaveMode{T,Acoustic{T,2},RadialSymmetry{2}}) where T

    # Unpacking parameters
    k = wavemode.ω / wavemode.source.medium.c
    k_eff = wavemode.wavenumber
    R = outer_radius(wavemode.material.shape)

    species = wavemode.material.microstructure.species
    rs = outer_radius.(species)

    N = wavemode.basis_order

    # Formula for the integral of [J_n(kr) * J_n(k_eff r) rdr] over the interval (0,R).
    I(n,R) = R/(k^T(2)-k_eff^T(2))*(k*besselj(n+1,k*R)*besselj(n,k_eff*R)-k_eff*besselj(n,k*R)*besselj(n+1,k_eff*R))

    # and an integral of the eigenvectors over the species
    Fscat0 = T(2pi)*sum(
        I(i[1]-1-N,R-rs[i[2]]) *
        wavemode.eigenvectors[i] *
        number_density(species[i[2]])
    for i in CartesianIndices(wavemode.eigenvectors))

    return [Fscat0]
end