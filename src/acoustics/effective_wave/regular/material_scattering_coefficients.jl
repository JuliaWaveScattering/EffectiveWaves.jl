kernelM(l::Int,x, y) = x * diffsbesselj(l,x) * sbesselj(l,y) - y * sbesselj(l,x) * diffsbesselj(l,y)

function material_scattering_coefficients(wavemode::EffectiveRegularWaveMode{T,3,Acoustic{T,3},WithoutSymmetry{3}}) where T

    # Unpacking parameters
    k = wavemode.ω / wavemode.medium.c

    species = wavemode.material.species
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

function material_scattering_coefficients(wavemode::EffectiveRegularWaveMode{T,3,Acoustic{T,3},AzimuthalSymmetry{3}}) where T

    # Unpacking parameters
    k = wavemode.ω / wavemode.medium.c

    species = wavemode.material.species
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
