using EffectiveWaves, Test, LinearAlgebra

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# ms = MultipleScattering # just in case Circle and Sphere conflict with definitions from other packages.

s1 = Specie(
    Acoustic(spatial_dim; ρ=10.2, c=10.1), Sphere(0.004);
    # Acoustic(spatial_dim; ρ=10.2, c=10.1), Sphere(0.4);
    volume_fraction=0.2
    # volume_fraction=0.01
);
s2 = Specie(
    Acoustic(spatial_dim; ρ=5.2, c=4.1), Sphere(0.003);
    # Acoustic(spatial_dim; ρ=5.2, c=4.1), Sphere(0.3);
    volume_fraction=0.1
);
species = [s1,s2]
species = [s1]

ω = 1e-1
# ω = 0.2
tol = 1e-7
basis_order = 1
# basis_order = 2

θ = 0.0
psource = PlaneSource(medium, [sin(θ),0.0,cos(θ)]);
source = plane_source(medium; direction = [sin(θ),0.0,cos(θ)])
material = Material(Sphere(4.0),species);

# x = rand(3)
# field(source,x,ω) - field(psource,x,ω)

R = outer_radius(material.shape)
# basis_field_order = estimate_regular_basisorder(typeof(source.medium), R * k_eff )
# basis_field_order = 3
basis_field_order = 3

k = ω / medium.c
ko = k * medium.c / s1.particle.medium.c

k_phi = wavenumber_low_volumefraction(ω, medium, species;
    basis_order = basis_order)

# AP_kps0 =  wavenumbers_bisection_robust(ω, medium, species;
#     bisection_mesh_points = 30,
#     basis_order = basis_order,
#     num_wavenumbers = 3, tol = tol,
#     symmetry = PlanarAzimuthalSymmetry())


opts = Dict(
    :tol => tol, :num_wavenumbers => 6,
    :mesh_size => 2.0, :mesh_points => 20,
    :basis_order => basis_order, :basis_field_order => basis_field_order
);
AP_kps = wavenumbers(ω, medium, species; symmetry = PlanarAzimuthalSymmetry(), opts...)

# wavemodes = WaveModes(ω, source, material; opts...)

k_eff = AP_kps[1]
# k_eff = A_kps[1]
(imag(k_eff) < -tol) && (k_eff = - k_eff)

eff_medium = effective_medium(medium, species)
k_low = ω/eff_medium.c;
eff_medium.c
ω / k_eff
ω / k_phi

setupsymmetry(source, material)
setupsymmetry(psource, material)

k_eff = AP_kps[2]
k_eff = AP_kps[1]

A_wave = WaveMode(ω, k_eff, psource, material;
    basis_order = basis_order, basis_field_order = basis_field_order);

scat_azi = material_scattering_coefficients(A_wave);

R_wave = WaveMode(ω, k_eff, source, material;
    basis_order = basis_order, basis_field_order = basis_field_order);

scat_reg = material_scattering_coefficients(R_wave);

norm(scat_azi - scat_reg) / norm(scat_azi)

Linc = basis_field_order + basis_order

effective_sphere = Particle(eff_medium, material.shape)

medium3 = Acoustic(spatial_dim; ρ=0.2, c=0.1)
medium3 = eff_medium
effective_sphere = Particle(medium3, material.shape)

Tmat = MultipleScattering.t_matrix(effective_sphere, medium, ω, Linc)
t_diag = [Tmat[i,i] for i in axes(Tmat,1)]


# medium2 = Acoustic(2; ρ=0.2, c=0.1)
# effective_sphere = Particle(medium2, Circle(outer_radius(material.shape)))
# outer_medium2 = Acoustic(2; ρ=1.0, c=1.0)
# Tmat = MultipleScattering.t_matrix(effective_sphere, outer_medium2, ω, Linc)

p = effective_sphere
outer_medium = medium
basis_order = Linc

n_to_l = [l for l = 0:Linc for m = -l:l];

inds_azi = findall([
    m == 0
for l = 0:Linc for m = -l:l]);

source_coefficients = regular_spherical_coefficients(source)(Linc,zeros(3),ω);

t_diag = [Tmat[i,i] for i in axes(Tmat,1)];

norm(scat_azi - t_diag[n_to_l .+ 1] .* source_coefficients) / norm(scat_azi)


MA = eigensystem(ω, psource, material;
    basis_order = basis_order,
    basis_field_order = basis_field_order
)

MM = eigensystem(ω, source, material;
    basis_order = basis_order,
    basis_field_order = basis_field_order
)

det(MM(k_eff))
det(MM(k_low))
det(MM(k_phi))


det(MA(k_eff))
det(MA(k_low))
det(MA(k_phi))

#NOTE: MM(k_eff) ≈ MM_svd.U * diagm(0 => MM_svd.S) * MM_svd.Vt
MM_svd = svd(MM(k_eff))
MM_svd.S
inds = findall(MM_svd.S .< 1e-4)
# 9 -> 48
eigvectors = MM_svd.V[:,inds]

MM_svd = svd(MA(k_eff))
MM_svd.S
inds = findall(MM_svd.S .< 1e-4)
# 3 -> 10
eigvectors = MM_svd.V[:,inds]

eigvectors = eigenvectors(ω, k_eff, source, material;
        basis_order = basis_order,
        basis_field_order = basis_field_order
)

eigvectors = eigenvectors(ω, k_eff, psource, material;
        basis_order = basis_order,
        basis_field_order = basis_field_order
)

sum(eigvectors[i] * α[i[2]] for i in CartesianIndices(eigvectors), dims = 2)




T = Float64

norm(MM(k_eff) * eigvectors[:,1] - MM_svd.S[inds[1]] .* eigvectors[:,1]) / norm(eigvectors[:,1])

norm(MM(k_eff) * eigvectors[:,end] - MM_svd.S[inds[end]] .* eigvectors[:,end]) / norm(eigvectors[:,end])

L = basis_order
L1 = basis_field_order
L2 = L1
# or is it this:
# L2 = L + L1

# l,dl <= L
# l1 <= L1
# l2 <= L1 + L
# l3 <= 2L # assuming L<= L1

# k1 = 1e-8
# k2 = 1e-8 + 1e-6im
# kernelN3D.(0:5, k1, k2) ./ (k2^2 - k1^2)
#
#
# outer_medium = medium
# p = s1.particle
# m = 0
#
# ak = outer_radius(p)*ω/outer_medium.c
#
# q = impedance(p.medium)/impedance(outer_medium) # Impedance ratio
# γ = outer_medium.c / p.medium.c #speed ratio
# numer = q * diffsbesselj(m, ak) * sbesselj(m, γ * ak) - sbesselj(m, ak)*diffsbesselj(m, γ * ak)
# denom = q * diffshankelh1(m, ak) * sbesselj(m, γ * ak) - shankelh1(m, ak)*diffsbesselj(m, γ * ak)
#
# - numer / denom
# t_matrix(s1.particle, medium, ω, 0)
#
#
# -(SphericalBesselJ[0, z]/(2 z)) +
#  1/2 (SphericalBesselJ[-1, z] - SphericalBesselJ[1, z])
#
# m=0
# (- sbesselj(m,ak) + (ak) * (sbesselj(m-1,ak) - sbesselj(m+1,ak))) / (2 * (ak))
