# EffectiveWaves.jl Documentation

*A Julia package for waves travelling in heterogeneous materials.*

The focus is on calculating the ensemble averaged waves, i.e. the statistical moments, of the waves. At present, the packages calculates effective wavenumbers, wave transmission and wave reflection from random particulate materials, see [preprint](https://arxiv.org/abs/1712.05427) and [notes](theory/WavesInMultiSpecies.pdf) for the formulas, and average wave scattering from spheres and cylinders filled with particles, see [preprint](https://arxiv.org/abs/2010.00934). The package can deal with any volume fraction and frequency range.

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


## Manual

You can learn to use this package through [examples](examples/README.md) or through our manual, which starts with a [Quick introduction](@ref).

## Contents
```@contents
Depth = 1
```
