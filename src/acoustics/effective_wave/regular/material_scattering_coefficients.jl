kernelM(l::Int,x, y) = x * diffsbesselj(l,x) * sbesselj(l,y) - y * sbesselj(l,x) * diffsbesselj(l,y)

function scattering_field(wavemode::EffectiveRegularWaveMode{T,3,Acoustic{T,3},WithoutSymmetry{3}}) where T
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

    v = regular_basis_function(wavemode.wavenumber, wavemode.medium)

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

function scattering_field(wavemode::EffectiveRegularWaveMode{T,3,Acoustic{T,3},Sym}) where {T, Sym<:AbstractAzimuthalSymmetry{3}}
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

    v = regular_basis_function(wavemode.wavenumber, wavemode.medium)

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
    scattering_field(wavemode::EffectiveRegularWaveMode{T,3,Acoustic{T,3},RadialSymmetry{3}})

In general ``F_{nn'} = \\delta_{m',-m} \\delta_{\\ell,\\ell'}(-1)^m F_{(\\ell,0),(\\ell,0)}``, so the scattering field is given by

``\\langle f_n \\rangle (x) = \\sum_{n'} F_{nn'} \\mathrm v_{n'}(k_* x) =  F_{(\\ell,0),(\\ell,0)} (-1)^m \\mathrm v_{(\\ell,-m)}(k_* x)``
"""
function scattering_field(wavemode::EffectiveRegularWaveMode{T,3,Acoustic{T,3},RadialSymmetry{3}}) where T
    L = wavemode.basis_order
    len = basisorder_to_basislength(Acoustic{T,3},L)

    ls, ms = spherical_harmonics_indices(L)
    inds = lm_to_spherical_harmonic_index.(ls,-ms)

    v = regular_basis_function(wavemode.wavenumber, wavemode.medium)

    function scat_field(x::AbstractVector{T})
        vs =  v(L, x)
        vs = vs[inds] .* (-one(T)) .^ ms

        vecs = [
            sum(wavemode.eigenvectors[ls[n]+1,i2,:]) * vs[n]
        for n = 1:len, i2 = 1:size(wavemode.eigenvectors,2)];

        return vecs
    end
end



function material_scattering_coefficients(wavemode::EffectiveRegularWaveMode{T,3,Acoustic{T,3},WithoutSymmetry{3}}) where T

    # Unpacking parameters
    k = wavemode.ω / wavemode.medium.c

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

function material_scattering_coefficients(wavemode::EffectiveRegularWaveMode{T,3,Acoustic{T,3},Sym}) where {T,Sym<:AbstractAzimuthalSymmetry{3}}

    # Unpacking parameters
    k = wavemode.ω / wavemode.medium.c

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

function material_scattering_coefficients(wavemode::EffectiveRegularWaveMode{T,3,Acoustic{T,3},RadialSymmetry{3}}) where T

    # Unpacking parameters
    k = wavemode.ω / wavemode.medium.c
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


function averaged_scattered_field(wavemodes::Array{EffectiveRegularWaveMode{T,2,Acoustic{T,2},RadialSymmetry{2}}}) where T

    # Unpacking parameters
    k = wavemodes[1].ω / wavemodes[1].medium.c
    R = outer_radius(wavemodes[1].material.shape)
    basis_order = wavemodes[1].basis_order

    # Compute the averaged coefficient of the material
    # only the one of zeroth order 'F0' is non zero

    # It requires to compute an integral of bessel functions Ipn for which we have a formula
    I(k_eff,n) = R/(k^2-k_eff)*(k*besselj(n+1,k*R)*besselj(n,k_eff*R)-k_eff*besselj(n,k*R)*besselj(n+1,k_eff*R)) # n ∈ [-basis_order:basis_order]

    # and an integral of the eigenvectors over the species
    F0 = 2pi*
    sum(
        sum(
            I(w.wavenumber,i[1]-1-basis_order) *
            w.eigenvectors[i] *
            number_density(w.material.microstructure.species[i[2]])
            for i in CartesianIndices(w.eigenvectors)
            )
        for w in wavemodes
        )


    return (x,y)-> F0*hankelh1(0,k*sqrt(x^2+y^2))
end
