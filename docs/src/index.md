# EffectiveWaves.jl Documentation

*A Julia package for calculating, processing and plotting waves travelling in heterogeneous materials. The focus is on ensemble averaged waves.*


At present, the package focuses on materails filled with randomly placed particles. You can calculate effective wavenumbers for 2D [[1](https://arxiv.org/abs/1712.05427)] and 3D [[4](https://arxiv.org/abs/2010.00934)] acoustics, wave transimission and wave reflection in 2D [[1](https://arxiv.org/abs/1712.05427),[2](https://arxiv.org/abs/1810.10816),[3](https://arxiv.org/abs/1905.06996)] and 3D [[4](https://arxiv.org/abs/2010.00934)], and scattering from an inhomogenious sphere [[4](https://arxiv.org/abs/2010.00934)]. See [these notes](theory/WavesInMultiSpecies.pdf) for brief formulas on effective wavenumbers.

Together with [MultipleScattering.jl](https://github.com/JuliaWaveScattering/MultipleScattering.jl), this package has been setup to easily extend to other dimensions, materials, and types of waves, such as elastic and electromagnetic waves.

!!! note
    First install Julia v1.0.5 or later, then run in the Julia REPL:
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

## References
[[1]](http://rspa.royalsocietypublishing.org/content/474/2212/20170864) Gower AL, Smith MJ, Parnell WJ, Abrahams ID. Reflection from a multi-species material and its transmitted effective wavenumber. Proc. R. Soc. A. 2018 Apr 1;474(2212):20170864.

[[2]](https://epubs.siam.org/doi/abs/10.1137/18M122306X) Gower, Artur L., William J. Parnell, and I. David Abrahams. "Multiple waves propagate in random particulate materials." SIAM Journal on Applied Mathematics 79.6 (2019): 2569-2592.

[[3]](https://royalsocietypublishing.org/doi/full/10.1098/rspa.2019.0344) Gower, Artur L., I. David Abrahams, and William J. Parnell. "A proof that multiple waves propagate in ensemble-averaged particulate materials." Proceedings of the Royal Society A 475.2229 (2019): 20190344.

[[4]](https://arxiv.org/abs/2010.00934) Gower, Artur Lewis, and Gerhard Kristensson. "Effective Waves for Random Three-dimensional Particulate Materials." arXiv preprint arXiv:2010.00934 (2020).

## Contents
```@contents
Depth = 1
```
