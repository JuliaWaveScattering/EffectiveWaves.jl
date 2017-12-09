# Multi-species effective waves

A Julia library for calculating, processing and plotting effective waves travelling in inhomogeneous materials.

At present, the library focuses on effective wavenumbers and wave reflection from random particulate materials, see ?? for details on the mathematics.

## Get started
This package is tested and works for Julia 0.6 and 0.5.
To get started, download and include the library
```julia
Pkg.clone("https://github.com/arturgower/EffectiveWaves.jl.git")
Pkg.add("SpecialFunctions")
```

## Simple example, complete code in [examples/demo.jl](examples/demo.jl)
### Run
Choose two types of particles, the first centred at [-2.,2.] and the second at [-2.,-2.]
```julia
include("src/multi-species.jl")

## Choose two species randomly (uniformly) distributed
# Usage Specie(ρ = density, r = radius, c = wavespeed, volfrac = volume fraction)
species = [
    Specie(ρ=WaterDistilled.ρ,r=30.e-6, c=WaterDistilled.c, volfrac=0.1),
    Specie(ρ=Inf, r=100.0e-6, c=2.0, volfrac=0.2)
]
# background medium
background = Glycerol
```
### A list of possible materials are given in [src/materials.jl](src/materials.jl).

### Calculate effective wavenumbers
```julia

# angular frequencies
ωs = linspace(0.01,1.0,60)*30.0e6
wavenumbers = sqrt.(multispecies_wavenumber(ωs, background, species))

speeds = ωs./real(wavenumbers)
attenuations = imag(wavenumbers)
```
