using EffectiveWaves, Test, LinearAlgebra

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# ms = MultipleScattering # just in case Circle and Sphere conflict with definitions from other packages.

s1 = Specie(
    Acoustic(spatial_dim; ρ=10.2, c=10.1), Sphere(0.4);
    volume_fraction=0.1
);
s2 = Specie(
    Acoustic(spatial_dim; ρ=5.2, c=4.1), Sphere(0.3);
    volume_fraction=0.1
);
species = [s1,s2]

ω = 0.5
tol = 1e-6
basis_order = 2

AP_det = dispersion_equation(ω, medium, species, PlanarAzimuthalSymmetry(); basis_order = basis_order)

k_eff = wavenumber_low_volumefraction(ω, medium, species;
    basis_order = basis_order)

AP_det(k_eff)

# AP_kps0 =  wavenumbers_bisection_robust(ω, medium, species;
#     bisection_mesh_points = 30,
#     basis_order = basis_order,
#     num_wavenumbers = 3, tol = tol,
#     symmetry = PlanarAzimuthalSymmetry())

AP_kps = wavenumbers(ω, medium, species;
    basis_order = basis_order,
    num_wavenumbers = 2, tol = tol,
    symmetry = PlanarAzimuthalSymmetry())

θ = 0.0
normal = [0.0,0.0,-1.0]; # an outward normal to the surface
material = Material(Sphere(4.0),species);
source2 = PlaneSource(medium, [sin(θ),0.0,cos(θ)]);
source = plane_source(medium; direction = [sin(θ),0.0,cos(θ)])

x = rand(3)
field(source,x,ω) - field(source2,x,ω)

basis_order = 2
basis_field_order = 4

k_eff = AP_kps[1]

setupsymmetry(source, material)

MM = eigensystem(ω, source, material;
    basis_order = basis_order,
    basis_field_order = basis_field_order
)

MM_svd = svd(MM(k_eff))

inds = findall(MM_svd.S .< tol)

L = basis_order
L1 = basis_field_order
L2 = L1
# or is it this:
# L2 = L + L1

# l,dl <= L
# l1 <= L1
# l2 <= L1 + L
# l3 <= 2L # assuming L<= L1



Freg = [
        (dl,dm,l1,m1)
    for dl = 0:L for dm = -dl:dl
for l1 = 0:L1 for m1 = -l1:l1]

reshape(Freg,(:,(L+1)^2))

Fazi = [
        (dl,dm,l1)
for dl = 0:L for dm = -dl:dl for l1 = abs(dm):L1]

reshape(Fazi,(:,(L+1)^2))
