using EffectiveWaves, Test, LinearAlgebra

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# ms = MultipleScattering # just in case Circle and Sphere conflict with definitions from other packages.

s1 = Specie(
    Acoustic(spatial_dim; ρ=10.2, c=10.1), Sphere(0.4);
    volume_fraction=0.2
    # volume_fraction=0.01
);
s2 = Specie(
    Acoustic(spatial_dim; ρ=5.2, c=4.1), Sphere(0.3);
    volume_fraction=0.1
);
species = [s1,s2]
species = [s1]

ω = 1e-5
# ω = 0.4
tol = 1e-7
basis_order = 1

k = ω / medium.c
ko = k * medium.c / s1.particle.medium.c

k_phi = wavenumber_low_volumefraction(ω, medium, species;
    basis_order = basis_order)

# AP_kps0 =  wavenumbers_bisection_robust(ω, medium, species;
#     bisection_mesh_points = 30,
#     basis_order = basis_order,
#     num_wavenumbers = 3, tol = tol,
#     symmetry = PlanarAzimuthalSymmetry())

AP_kps = wavenumbers(ω, medium, species;
    basis_order = basis_order,
    num_wavenumbers = 2, tol = tol,
    mesh_size = 2.0,
    mesh_points = 10,
    symmetry = PlanarAzimuthalSymmetry())

θ = 0.0
material = Material(Sphere(4.0),species);
psource = PlaneSource(medium, [sin(θ),0.0,cos(θ)]);
source = plane_source(medium; direction = [sin(θ),0.0,cos(θ)])

x = rand(3)
field(source,x,ω) - field(psource,x,ω)

k_eff = AP_kps[1]
(imag(k_eff) < -tol) && (k_eff = - k_eff)

ρ = medium.ρ
β = medium.c^2 * medium.ρ

# φ = sum(volume_fraction.(species));
φ = volume_fraction(species[1]);
ρo = s1.particle.medium.ρ
βo = s1.particle.medium.c^2 * s1.particle.medium.ρ

Dρ = (ρ - ρo)/(ρ + 2.0*ρo)
Dβ = (β - βo)/βo
βp = ( 1.0/β + Dβ/β * φ )^(-1)
ρp = ρ * (1.0 - Dρ * φ)/(1.0 + 2.0 * Dρ * φ)
      # ρ*(1 - ρ_frac)/(1 + T(Dim - 1) * ρ_frac)

cp = sqrt(
    (1/β + Dβ/β * φ)^(-1) * (ρ * (1 - Dρ * φ)/( 1 + 2.0 * Dρ * φ ))^(-1)
)

eff_medium = effective_medium(medium, species)
k_low = ω/eff_medium.c;
eff_medium.c
ω / k_eff
ω / k_phi

R = outer_radius(material.shape)
basis_field_order = estimate_regular_basisorder(typeof(source.medium), R * k_eff )
basis_field_order = 2

setupsymmetry(source, material)
setupsymmetry(psource, material)

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
eigvectors = MM_svd.V[:,inds]

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

k1 = 1e-8
k2 = 1e-8 + 1e-6im
kernelN3D.(0:5, k1, k2) ./ (k2^2 - k1^2)


outer_medium = medium
p = s1.particle
m = 0

ak = outer_radius(p)*ω/outer_medium.c

q = impedance(p.medium)/impedance(outer_medium) # Impedance ratio
γ = outer_medium.c / p.medium.c #speed ratio
numer = q * diffsbesselj(m, ak) * sbesselj(m, γ * ak) - sbesselj(m, ak)*diffsbesselj(m, γ * ak)
denom = q * diffshankelh1(m, ak) * sbesselj(m, γ * ak) - shankelh1(m, ak)*diffsbesselj(m, γ * ak)

- numer / denom
t_matrix(s1.particle, medium, ω, 0)


-(SphericalBesselJ[0, z]/(2 z)) +
 1/2 (SphericalBesselJ[-1, z] - SphericalBesselJ[1, z])

m=0
(- sbesselj(m,ak) + (ak) * (sbesselj(m-1,ak) - sbesselj(m+1,ak))) / (2 * (ak))
