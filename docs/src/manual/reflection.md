# Reflection and transmission

Once we have calculated the effective wavenumber, we can calculate the average reflection and transmission.

Here we will use an incident plane wave

$u_\text{in} = U e^{i \mathbf k \cdot \mathbf r},$

where $U$ is the amplitude, which is usually set to $U=1$ for acoustics, $k = \|\mathbf k\|$ is the incident wavenumber, and $\mathbf k$ can be a two or three dimensional vector.

We will assume that $u_\text{in}$ is arriving from $z<0$ ($x<0$ for 2D). If the material occupies the region $\mathcal R = \{ z>0 : \mathbf r \in \mathbb R^3\}$, then the reflected wave will be given by  

$u_\text{R} = R e^{i (k_x x + k_y y - k_z z)}.$

The code below calculates $R$, which is called the reflection coefficient. Both reflection and transmission are simpler to calculate when there exists only [one effective wavenumber](@ref two-dim-acoustics-one_reflection). Currently, we have only implemented the reflection coefficient for multiple effective wavenumbers for [2D acoustics](@ref two-dim-acoustic-multiple-reflection).

# [2D acoustics](@id two-dim-acoustics-reflection)

## [Low frequency reflection](@id two-dim-acoustic-low-reflection)

The simplest case is for low frequency, where the average reflection coefficient $R$ reduces to the reflection coefficient from a homogeneous material:

$R = \frac{q_* \cos \theta_\text{in} - \cos \theta_*}{q_* \cos \theta_\text{in} + \cos \theta_*},$

where

$\mathcal R = \{ x>0 : \mathbf r \in \mathbb R^2\}  \;\; \text{and} \;\; \mathbf k = k (\cos \theta_\text{in}, \sin \theta_\text{in}).$

In code this becomes

```julia 2; setup = :(using EffectiveWaves)
spatial_dim = 2
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# Choose the species
radius1 = 0.001
s1 = Specie(
    Acoustic(spatial_dim; ρ=10.2, c=10.1), radius1;
    volume_fraction=0.2
);
radius2 = 0.002
s2 = Specie(
    Acoustic(spatial_dim; ρ=0.2, c=4.1), radius2;
    volume_fraction=0.15

);
species = [s1,s2]

# Choose the frequency
ω = 1e-2
k = ω / real(medium.c)

# Check that we are indeed in a low frequency limit
k * radius2 < 1e-4

eff_medium = effective_medium(medium, species)

normal = [-1.0,0.0] # an outward normal to the surface

# Define the material region
material = Material(Halfspace(normal),species)

# define a plane wave source travelling at a 45 degree angle in relation to the material
source = PlaneSource(medium, [cos(pi/4.0),sin(pi/4.0)])

R = reflection_coefficient(ω,source, eff_medium, material.shape)

# output

0.13666931757028777 + 0.0im
```

## [One plane wave mode](@id two-dim-acoustics-one_reflection)

Note that for there are formulas for low volume fraction expansions of the reflection coefficient, see [`reflection_coefficient_low_volumefraction`](@ref). However, it is almost better to use the exact expression, as the added computational cost is minimal.  

As an example, we will use the same material defined for the [low frequency case](@ref two-dim-acoustic-low-reflection)
```julia 2

k_effs = wavenumbers(ω, medium, species; tol = 1e-6, num_wavenumbers = 1, basis_order = 1)

# Calculate the wavemode for the first wavenumber
wave1 = WaveMode(ω, k_effs[1], source, material; tol = 1e-6, basis_order = 1)

R = reflection_coefficient(ω, wave1, source, material)

# output

0.13666931756494047 - 5.127394485569188e-14im
```

## [Multiple effective modes](@id two-dim-acoustic-multiple-reflection)

See [A numerical matching method](@ref) for an example that uses multiple effective wave modes to calculate the reflection coefficient.

# 3D acoustics

## Low frequency reflection

The simplest case is for low frequency, where the average reflection coefficient $R$ reduces to the reflection coefficient from a homogeneous material:

$R = \frac{q_* \cos \theta_\text{in} - \cos \theta_*}{q_* \cos \theta_\text{in} + \cos \theta_*},$

where

$\mathcal R = \{ z>0 : \mathbf r \in \mathbb R^3\} \;\; \text{and} \;\; \mathbf k = k (\cos \phi_\text{in} \sin \theta_\text{in}, \sin \phi_\text{in} \sin \theta_\text{in}, \cos \theta_\text{in}).$

In code this becomes
```julia 2; setup = :(using EffectiveWaves)

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# Choose the species
radius1 = 0.001
s1 = Specie(
    Acoustic(spatial_dim; ρ=10.2, c=10.1), radius1;
    volume_fraction=0.2
);
radius2 = 0.002
s2 = Specie(
    Acoustic(spatial_dim; ρ=0.2, c=4.1), radius2;
    volume_fraction=0.15

);
species = [s1,s2]

# Choose the frequency
ω = 1e-2
k = ω / real(medium.c)

# Check that we are indeed in a low frequency limit
k * radius2 < 1e-4

eff_medium = effective_medium(medium, species)

normal = [0.0,0.0,-1.0] # an outward normal to the surface

# Define the material region
material = Material(Halfspace(normal),species)

# define a plane wave source travelling at a 45 degree angle in relation to the material
source = PlaneSource(medium, [cos(pi/4.0),0.0,sin(pi/4.0)])

R = reflection_coefficient(ω,source, eff_medium, material.shape)

# output

0.15901291515072696 + 1.9230884702465903e-17im
```

## [One plane wave mode ](@id three-dim-acoustic-one-reflection)

Using only one plane wave mode we can calculate both reflection and transmission from a plate.

```julia 2
k_effs = wavenumbers(ω, medium, species; tol = 1e-6, num_wavenumbers = 1, basis_order = 1)

# Define a plate
normal = [0.0,0.0,-1.0] # an outward normal to both surfaces of the plate
width = 1.0 # plate width

# Define the material region
material = Material(Plate(normal,width),species)

basis_order = 1;
kws = Dict(:basis_order => basis_order)

MM = eigensystem(ω, source, material; kws...)

# calculate eigenvectors
using LinearAlgebra

k_eff = k_effs[1]
MM_svd = svd(MM(k_eff))
MM_svd.S[end]
v1 = MM_svd.V[:,end]

k_eff = -k_effs[1]
MM_svd = svd(MM(k_eff))
MM_svd.S[end]
v2 = MM_svd.V[:,end]




vecs = eigenvectors(ω, k_effs[1], source, material; kws...)


# Calculate the wavemode for the first wavenumber
wave1 = WaveMode(ω, k_effs[1], source, material; tol = 1e-6, basis_order = 1)


```
Currently implementing... formulas from [Gower & Kristensson 2020](https://arxiv.org/pdf/2010.00934.pdf).
