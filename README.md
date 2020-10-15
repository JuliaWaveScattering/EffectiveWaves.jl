# EffectiveWaves

*A Julia package for calculating, processing and plotting waves travelling in heterogeneous materials. The focus is on ensemble averaged waves.*

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
|[![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url] [![][codecov-img]][codecov-url] [![][coveralls-img]][coveralls-url] |


<!-- You can run Julia on [JuliaBox](https://www.juliabox.com/) in your browser without installation. -->

At present, the package focuses on materails filled with randomly placed particles. You can calculate effective wavenumbers for 2D [1](https://arxiv.org/abs/1712.05427) and 3D [2](https://arxiv.org/abs/2010.00934) acoustics, wave transimission and wave reflection in 2D [1](https://arxiv.org/abs/1712.05427)[3](https://arxiv.org/abs/1810.10816)[4](https://arxiv.org/abs/1905.06996) and 3D [2](https://arxiv.org/abs/2010.00934), and scattering from an inhomogenious sphere [2](https://arxiv.org/abs/2010.00934). See [these notes](docs/src/theory/WavesInMultiSpecies.pdf) for brief formulas on effective wavenumbers.


## Installation
Type into Julia:
```julia
using Pkg
Pkg.clone("https://github.com/arturgower/EffectiveWaves.jl.git")

using EffectiveWaves
```

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **documentation of the most recently tagged version.**
- [**DEVEL**][docs-dev-url] &mdash; *documentation of the in-development version.*

## Simple example
Effective wavenumbers for two species randomly (uniformly) distributed in Glycerol.
```julia
#where: ρ = density, r = radius, c = wavespeed, and volfrac = volume fraction

const WaterDistilled= Medium(ρ=0.998*1000, c = 1496.0)
const Glycerol      = Medium(ρ=1.26*1000,  c = 1904.0)

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
wavenumbers = wavenumber_low_volumefraction(ωs, background, species)

speeds = ωs./real(wavenumbers)
attenuations = imag(wavenumbers)
```
For a list of possible materials go to [src/materials.jl](src/materials.jl).

## More examples
For more examples and details go to [docs/src/examples/](docs/src/examples/).

## Acknowledgements and contributing
This library was originally written by [Artur L Gower](https://arturgower.github.io/).
Contributions are very welcome.

## References
[[1]](http://rspa.royalsocietypublishing.org/content/474/2212/20170864) Gower AL, Smith MJ, Parnell WJ, Abrahams ID. Reflection from a multi-species material and its transmitted effective wavenumber. Proc. R. Soc. A. 2018 Apr 1;474(2212):20170864.

[[2]](https://epubs.siam.org/doi/abs/10.1137/18M122306X) Gower, Artur L., William J. Parnell, and I. David Abrahams. "Multiple waves propagate in random particulate materials." SIAM Journal on Applied Mathematics 79.6 (2019): 2569-2592.

[[3]](https://royalsocietypublishing.org/doi/full/10.1098/rspa.2019.0344) Gower, Artur L., I. David Abrahams, and William J. Parnell. "A proof that multiple waves propagate in ensemble-averaged particulate materials." Proceedings of the Royal Society A 475.2229 (2019): 20190344.


[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JuliaWaveScattering.github.io/EffectiveWaves.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaWaveScattering.github.io/EffectiveWaves.jl/stable

[travis-img]: https://travis-ci.org/JuliaWaveScattering/EffectiveWaves.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaWaveScattering/EffectiveWaves.jl

[codecov-img]: http://codecov.io/github/JuliaWaveScattering/EffectiveWaves.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaWaveScattering/EffectiveWaves.jl?branch=master

[coveralls-img]: https://coveralls.io/repos/github/JuliaWaveScattering/EffectiveWaves.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/JuliaWaveScattering/EffectiveWaves.jl?branch=master

[issues-url]: https://github.com/JuliaWaveScattering/EffectiveWaves.jl/issues
