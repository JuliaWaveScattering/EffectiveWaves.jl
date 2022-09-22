# Reflection and transmission

Once we have calculated the effective wavenumber, we can calculate the average reflection and transmission.

Here we will use an incident plane wave

$u_\text{in} = U e^{i \mathbf k \cdot \mathbf r},$

where $U$ is the amplitude, which is usually set to $U=1$ for acoustics, $k = \|\mathbf k\|$ is the incident wavenumber, and $\mathbf k$ can be a two or three dimensional vector.

We will assume that $u_\text{in}$ is arriving from $z<0$ ($x<0$ for 2D). If the material occupies the region $\mathcal R = \{ z>0 : \mathbf r \in \mathbb R^3\}$, then the reflected wave will be given by  

$u_\text{R} = R e^{i (k_x x + k_y y - k_z z)}.$

The code below calculates $R$, which is called the reflection coefficient. Both reflection and transmission are simpler to calculate when there exists only [one effective wavenumber](@ref two-dim-acoustics-one_reflection). Currently, we have only implemented the reflection coefficient for multiple effective wavenumbers for [2D acoustics](@ref two-dim-acoustic-multiple-reflection).

Many tests for 3D reflection and transmission are in [test/effective-3D/planar-symmetry.jl](../../../test/effective-3D/planar-symmetry.jl), while tests for 2D are in the files in the folder [test/effective-2D](../../../test/effective-2D).  

The formulas used to calculate the below can mostly be found in [Gower & Kristensson 2020](https://arxiv.org/pdf/2010.00934.pdf) and [Gower et al. 2018](https://arxiv.org/abs/1712.05427)

# [2D acoustics](@id two-dim-acoustics-reflection)

## [Low frequency reflection](@id two-dim-acoustic-low-reflection)

The simplest case is for low frequency, where the average reflection coefficient $R$ reduces to the reflection coefficient from a homogeneous material:

$R = \frac{q_* \cos \theta_\text{in} - \cos \theta_*}{q_* \cos \theta_\text{in} + \cos \theta_*},$

where

$\mathcal R = \{ x>0 : \mathbf r \in \mathbb R^2\}  \;\; \text{and} \;\; \mathbf k = k (\cos \theta_\text{in}, \sin \theta_\text{in}).$

In code this becomes

```jldoctest 2; setup = :(using EffectiveWaves)
spatial_dim = 2
medium = Acoustic(spatial_dim; ρ=1.2, c=1.5)

# Choose the species
radius1 = 0.1
s1 = Specie(
    Acoustic(spatial_dim; ρ=10.2, c=10.1), radius1;
    volume_fraction=0.2
);
radius2 = 0.2
s2 = Specie(
    Acoustic(spatial_dim; ρ=0.2, c=4.1), radius2;
    volume_fraction=0.15

);
species = [s1,s2]

# Choose the frequency
ω = 1e-4
k = ω / medium.c

# Calculate the equivalent effective medium in the asymptotic low frequency limit
eff_medium = effective_medium(medium, species)

normal = [-1.0,0.0] # an outward normal to the surface

# Define the material region
material = Material(Halfspace(normal),species)

# define a plane wave source travelling at a 45 degree angle in relation to the material
source = PlaneSource(medium, [cos(pi/4.0),sin(pi/4.0)])

R = reflection_coefficient(ω,source, eff_medium, material.shape);

round(R * 100) / 100
# output

0.14 + 0.0im
```

## [One plane wave mode](@id two-dim-acoustics-one_reflection)

Note there are formulas for low volume fraction expansions of the reflection coefficient, see [`reflection_coefficient_low_volumefraction`](@ref). However, it is better to use the exact expression, as the different in computational cost is minimal.  

As an example, we will use the same material defined for the [low frequency case](@ref two-dim-acoustic-low-reflection)
```jldoctest 2

k_effs = wavenumbers(ω, medium, species; tol = 1e-6, num_wavenumbers = 1, basis_order = 1)

# Calculate the wavemode for the first wavenumber
wave1 = WaveMode(ω, k_effs[1], source, material; tol = 1e-6, basis_order = 1)

R = reflection_coefficient(ω, wave1, source, material)

round(R * 100) / 100

# output
0.14 + 0.0im
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
```jldoctest 2

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=0.3, c=0.5)

# Choose the species
radius1 = 0.1
s1 = Specie(
    Acoustic(spatial_dim; ρ=10.2, c=10.1), radius1;
    volume_fraction=0.2
);
species = [s1]

# Choose the frequency
ω = 1e-5
k = ω / medium.c

# For the limit of low frequencies we can define
eff_medium = effective_medium(medium, species)

# Define a plate
r = maximum(outer_radius.(species))
normal = [0.0,0.0,-1.0] # an outward normal to both surfaces of the
width = 150.0 # plate width
origin = [0.0,0.0,width/2] # the centre of the plate


# the size of the effective low frequency limit material is one particle radius smaller
plate_low = Plate(normal,width - 2r,origin)
halfspace_low = Halfspace(normal,halfspace.origin - r)

# define a plane wave source travelling at a 45 degree angle in relation to the material
source = PlaneSource(medium, [cos(pi/4.0),0.0,sin(pi/4.0)])

Ramp1 = reflection_coefficient(ω, source, eff_medium, halfspace_low)

# planewave_amplitudes returns the
amps = planewave_amplitudes(ω, source, eff_medium, plate_low)
Ramp = amps[1]
Tamp = amps[2]

round(1000 * Tamp) / 1000

# output

1.0 + 0.002im
```
The function [`planewave_amplitudes`](@ref) returns `[R, T, P1, P2]` where `R` is the reflection coefficient, `T` is the coefficient of the transmitted wave, and `P1` (`P2`) are the amplitudes of the wave travling forward (backward) inside the plate.


## [One plane wave mode ](@id three-dim-acoustic-one-reflection)

Using only one plane wave mode we can calculate both reflection and transmission from a plate.

```julia 2
k_effs = wavenumbers(ω, medium, species;
    tol = 1e-6,
    num_wavenumbers = 1,
    basis_order = 1
)
k_eff = k_effs[1]

abs(k_eff - ω / eff_medium.c) < 1e-10

halfspace = Halfspace(normal)
plate = Plate(normal,width,origin)

material = Material(halfspace,species)

# Calculate the wavemode for the first wavenumber
# the WaveMode function calculates the types of waves and solves the needed boundary conditions
wavemode = WaveMode(ω, k_eff, source, material; tol = 1e-6, basis_order = 1);

Reff = reflection_coefficient(wavemode, source, material)

material = Material(plate,species)
wavemodes = WaveMode(ω, k_eff, source, material; tol = 1e-6, basis_order = 1);
RTeff = reflection_transmission_coefficients(wavemodes, source, material);

abs(Ramp1 - Reff) < 1e-6

abs.(RTeff - [Ramp; Tamp]) .< [1e-4, 5e-4]

# Note that summing the reflection and transmission from the homogeneous low frequency medium gives 1
sum(abs.([Ramp,Tamp]).^2) ≈ 1.0
sum(abs.(RTeff).^2)
```
Formulas from [Gower & Kristensson 2020](https://arxiv.org/pdf/2010.00934.pdf).
