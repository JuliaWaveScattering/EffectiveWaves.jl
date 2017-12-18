# Multi-species effective waves

A Julia (v0.5-v0.6) library for calculating, processing and plotting effective waves travelling in inhomogeneous materials.
You can run Julia on [JuliaBox](https://www.juliabox.com/) in your browser without installation.

At present, the library focuses on effective wavenumbers and wave reflection from random particulate materials, see [arXiv preprint](https://arxiv.org/abs/1712.05427) for details on the mathematics, or [these notes](theory/MultispeciesWaves.pdf) for the formulas.

## Get started
Type into Julia:
```julia
Pkg.clone("https://github.com/arturgower/EffectiveWaves.jl.git")
using EffectiveWaves
```

## Simple example
Effective wavenumbers for two species randomly (uniformly) distributed in Glycerol, complete code in [examples/demo.jl](examples/demo.jl).
```julia
#where: ρ = density, r = radius, c = wavespeed, and volfrac = volume fraction

species = [
    Specie(ρ=WaterDistilled.ρ,r=30.e-6, c=WaterDistilled.c, volfrac=0.1),
    Specie(ρ=Inf, r=100.0e-6, c=2.0, volfrac=0.2)
]
# background medium
background = Glycerol
```

Calculate effective wavenumbers:
```julia

# angular frequencies
ωs = linspace(0.01,1.0,60)*30.0e6
wavenumbers = multispecies_wavenumber(ωs, background, species)

speeds = ωs./real(wavenumbers)
attenuations = imag(wavenumbers)
```
For a list of possible materials go to [examples/materials.jl](examples/materials.jl).

## More examples
For more examples and details go to [examples/](examples/).

## Acknowledgements and contributing
This library was originally written by Artur L Gower.
Please contribute, if nothing else, criticism is welcome, as I am relatively new to Julia.
