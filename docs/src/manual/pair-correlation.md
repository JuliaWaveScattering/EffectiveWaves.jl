# Pair correlation

In the previous example there was no mention on how the particles are distributed. This can be specified by choosing a pair correlation. Choosing the particle distribution only affects the effective wavenumbers and wavemodes, see refs...

## Percus-Yevick

Let us consider a material filled with only one type of particle

```jldoctest pair; setup = :(using EffectiveWaves), output = false, filter = r".*"s
medium = Acoustic(3; ρ=1.2, c=1.0)

# Choose the species
r = 0.5
s = Specie(
    Acoustic(3; ρ = 0.1, c = 0.1),
    Sphere(r),
    volume_fraction = 0.3,
    exclusion_distance = 1.01
);

# output

```
Next we create a microstructure that has only this species, and has a specific pair-correlation

```jldoctest pair; output = false, filter = r".*"s

pair_type = PercusYevick(rtol = 1e-2, maxsize = 200)

micro = Microstructure(s, pair_type);

# output
```
and we can plot the result
```julia
using Plots

plot(micro.paircorrelations[1].r, 1.0 .+ micro.paircorrelations[1].dp,
    xlab = "distance", ylab = "P-Y"
)
```
![../PY-30-pair.png](../assets/PY-30-pair.png)

which we can compare with Figure 8.3.1 from [1] below.

![../TKD-PY-30.jpg](../assets/TKD-PY-30.jpg)

Note that for $x < 1$ the two particles of radius 0.5 would overlap, so the pair correlation should be zero. Also note that `dp` is the variation from uncorrelated, which is why we add 1.0 to get the pair correlation.

# Calculate an effective wavenumber

The more points sampled within the pair correlation the longer it will take to calculate the effective wavenumber.

First we calculate the wavenumbers with the simplest pair correlation (hole correction), and then compare the results with Percus-Yevick.

```jldoctest pair; output = false, filter = r".*"s

micro = Microstructure(s);

ω = 2.0

kps = wavenumbers(ω, medium, micro;
    basis_order = 1, num_wavenumbers = 5
)

pair_type = PercusYevick(rtol = 5e-2, maxsize=40)
micro = Microstructure(s, pair_type);

kps2 = wavenumbers(ω, medium, micro;
    basis_order = 1, num_wavenumbers = 5
)

```

## References

[1] Kong, Jin Au, Leung Tsang, Kung-Hau Ding, and Chi On Ao. Scattering of electromagnetic waves: numerical simulations. John Wiley & Sons, 2004.

[[2]](https://github.com/JuliaWaveScattering/EffectiveWaves.jl/blob/master/docs/src/theory/P-Y.pdf) Gerhard Kristensson. "The Percus-Yevick approximation". [github.com/JuliaWaveScattering/EffectiveWaves.jl][https://github.com/JuliaWaveScattering/EffectiveWaves.jl] (2022).
