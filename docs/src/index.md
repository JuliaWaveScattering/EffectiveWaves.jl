# EffectiveWaves.jl Documentation

A Julia package for calculating, processing and plotting waves travelling in heterogeneous materials. The focus is on calculating the ensemble averaged waves, i.e. the statistical moments, of the waves.
You can run Julia on [JuliaBox](https://www.juliabox.com/) in your browser without installation.

At present, the packages calculates effective wavenumbers, wave transimission and wave reflection from random particulate materials in two-dimensions, see [arXiv preprint](https://arxiv.org/abs/1712.05427) for some details on the mathematics, or [these notes](theory/MultispeciesWaves.pdf) for the formulas.

!!! note
    First install Julia v1.0.0 or later, then run in the Julia REPL:
    ```julia
    using Pkg
    Pkg.add EffectiveWaves
    ```
    Alternatively, use the Julia package manager.
    From the Julia REPL, type `]` to enter the Pkg REPL mode and then run
    ```
    pkg> add EffectiveWaves
    ```

## Quick introduction

Here we calculate the effective wavenumbers for two species randomly (uniformly) distributed in Glycerol, complete code in [examples/demo.jl](examples/demo.jl).
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
ωs = LinRange(0.01,1.0,60)*30.0e6
wavenumbers = wavenumber_low_volfrac(ωs, background, species)

speeds = ωs./real(wavenumbers)
attenuations = imag(wavenumbers)
```
For a list of possible materials go to [examples/materials.jl](examples/materials.jl).

## Contents
```@contents
Depth = 1
```

## Index

```@index
```
