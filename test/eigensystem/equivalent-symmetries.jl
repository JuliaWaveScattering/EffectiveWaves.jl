
@testset "symmetry equivalent wavenumbers and eigenvectors" begin
using EffectiveWaves, Test, LinearAlgebra

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# ms = MultipleScattering # just in case Circle and Sphere conflict with definitions from other packages.

s1 = Specie(
    Acoustic(spatial_dim; ρ=0.2, c=0.1), Sphere(0.4);
    volume_fraction=0.2
);
species = [s1]

ω = 0.9
tol = 1e-7

AP_det = dispersion_equation(ω, medium, species, PlanarAzimuthalSymmetry())
P_det = dispersion_equation(ω, medium, species, PlanarSymmetry())
AR_det = dispersion_equation(ω, medium, species, AzimuthalSymmetry())
R_det = dispersion_equation(ω, medium, species, WithoutSymmetry())

AP_kps = wavenumbers(ω, medium, species;
    num_wavenumbers = 6, tol = tol,
    symmetry = PlanarAzimuthalSymmetry())

P_kps = wavenumbers(ω, medium, species;
    num_wavenumbers = 6, tol = tol,
    symmetry = PlanarSymmetry())

# The dispersion equations R_det become very unstable when using effective wavenumbers with higher dispersion.
AP_kps = AP_kps[1:min(length(AP_kps),6)]
P_kps = P_kps[1:min(length(P_kps),6)]

# As plane waves with azimuthal symmetry is a sub-case of plane-waves, and all materials allow for the effective wavenumbers of plane waves, all the below determinant equations should be satisfied
@test maximum(AP_det.(AP_kps)) < tol
@test maximum(P_det.(AP_kps)) < tol
@test maximum(AR_det.(AP_kps)) < tol^2
@test maximum(R_det.(AP_kps)) < tol^3

P_disp = dispersion_complex(ω, medium, species,  PlanarSymmetry(); tol = tol)

# However, there do exist effective wavenumbers for plane-waves which have eigen-vectors that do not satisfy azimuthal symmetry. This is why maximum(AP_det.(P_kps)) != 0.0
@test maximum(P_det.(P_kps)) < tol
@test maximum(AR_det.(P_kps)) < tol^2
@test maximum(R_det.(P_kps)) < tol^3

# Choose specific materials to calculate eignvectors

basis_order = 2
basis_field_order = 2*basis_order

R_MM = eigensystem(ω, medium, species,WithoutSymmetry();
        basis_order = basis_order,
        basis_field_order = basis_field_order
)

k_eff = P_kps[1]
MM_svd = svd(R_MM(k_eff))
inds = findall(abs.(MM_svd.S) .< 1e-4)

Rvs = [MM_svd.V[:,i] for i in inds] # eigenvectors

θp = 0.2
φp = 0.1

Ys = spherical_harmonics(basis_field_order, θp, φp);
ls, ms = spherical_harmonics_indices(basis_field_order)


L = basis_order
L1 = basis_field_order

# data = [ [dl,dm,l1,m1] for dl = 0:L for dm = -dl:dl for l1 = 0:L1 for m1 = -l1:l1];
# data = transpose(reshape(data, :, (L+1)^2))
#
# [ [data[i] i[2]]  for i in CartesianIndices(data) ]

RvMs = [transpose(reshape(v, :, (basis_order+1)^2)) for v in Rvs]
# Rvs[1], :, (basis_order+1)^2))

Pvs = [
    sum([ vM[i] * Ys[i[2]] / 1.0im^ls[i[2]] for i in CartesianIndices(vM)], dims=2)[:]
for vM in RvMs]


P_MM = eigensystem(ω, medium, species, PlanarSymmetry();
        θp = θp, φp = φp,
        basis_order = basis_order
)

M = P_MM(k_eff)

@test maximum(norm(M * Pv) for Pv in Pvs) < 1e-7

end
