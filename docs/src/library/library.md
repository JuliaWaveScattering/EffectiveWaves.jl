# Base

```@meta
CurrentModule = EffectiveWaves
```

```@contents
Pages = ["base.md"]
```

## Defining the material

```@docs
PhysicalMedium
Acoustic
Specie
Material
Source
```
Currently the physical medium can be [Acoustic](https://juliawavescattering.github.io/MultipleScattering.jl/dev/library/acoustics/)s. Both Acoustic, PlaneSource, and Source are all imported from the package [MultipleScattering](https://juliawavescattering.github.io/MultipleScattering.jl/dev/). In the future these will be moved to a new package WaveScatteringBase.
Much of the code is dispatched based on the underlying symmetries of the problem

## The symmetry of the material and source
The symmetry shared between the material shape and source are used to specialise the form of the wavemode, see [Background](@ref).  
```@docs
setupsymmetry
WithoutSymmetry
PlanarSymmetry
AzimuthalSymmetry
PlanarAzimuthalSymmetry
```


## Types of waves

There are two main types used.

```@docs
EffectivePlaneWaveMode
EffectiveRegularWaveMode
```

## Effective wavenumbers and wavemodes

```@docs
wavenumbers
wavenumber_low_volumefraction
effective_medium
```

## Effective wavemodes and eignvectors

```@docs
WaveMode
WaveModes
```

## Reflection coefficients

```@docs
reflection_coefficient
reflection_coefficient(::T, ::EffectivePlaneWaveMode{T}, ::PlaneSource{T,2,1,Acoustic{T,2}}, ::Material{2,Halfspace{T,2}}) where T<:AbstractFloat
reflection_coefficient_low_volumefraction
```

## Transmission

```@docs
transmission_direction
transmission_angle(::Vector,::Vector)
transmission_angle(::SVector{2,CT} where CT <: Union{T,Complex{T}}, ::SVector{2,T}) where {T<:AbstractFloat}
transmission_angle(::SVector{3,CT} where CT <: Union{T,Complex{T}}, ::SVector{3,T}) where {T<:AbstractFloat}
```
