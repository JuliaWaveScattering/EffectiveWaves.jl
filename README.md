# EffectiveWaves

*A Julia package for calculating, processing and plotting waves travelling in heterogeneous materials. The focus is on ensemble averaged waves.*

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url] [![][codecov-img]][codecov-url] [![][coveralls-img]][coveralls-url] |


<!-- You can run Julia on [JuliaBox](https://www.juliabox.com/) in your browser without installation. -->

At present, the packages calculates effective wavenumbers, wave transimission and wave reflection from random particulate materials in two-dimensions, see [arXiv preprint](https://arxiv.org/abs/1712.05427) for details on the mathematics, or [these notes](docs/src/theory/WavesInMultiSpecies.pdf) for the formulas.

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
Please contribute, if nothing else, criticism is welcome, as I am relatively new to Julia.

## References
[[1]](http://rspa.royalsocietypublishing.org/content/474/2212/20170864) Gower AL, Smith MJ, Parnell WJ, Abrahams ID. Reflection from a multi-species material and its transmitted effective wavenumber. Proc. R. Soc. A. 2018 Apr 1;474(2212):20170864.

[[2]](https://arxiv.org/abs/1712.05427) Gower AL, Smith MJ, Parnell WJ, Abrahams ID. Reflection from a multi-species material and its transmitted effective wavenumber. arXiv preprint arXiv:1712.05427. 2017 Dec.

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://arturgower.github.io/EffectiveWaves.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://arturgower.github.io/EffectiveWaves.jl/stable

[travis-img]: https://travis-ci.org/arturgower/EffectiveWaves.jl.svg?branch=master
[travis-url]: https://travis-ci.org/arturgower/EffectiveWaves.jl

[codecov-img]: http://codecov.io/github/arturgower/EffectiveWaves.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/arturgower/EffectiveWaves.jl?branch=master

[coveralls-img]: https://coveralls.io/repos/github/arturgower/EffectiveWaves.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/arturgower/EffectiveWaves.jl?branch=master

[issues-url]: https://github.com/arturgower/EffectiveWaves.jl/issues
