# A numerical matching method
Both of these methods calculate the ensemble average wave in a particulate material. To do so, we need to specify the geometry of the material, as this ultimately determines the fields. We choose the simplest case which is a halfspace x>0 filled with particles.

## Calculate the effective wavemodes
```julia
using EffectiveWaves

## choose a 2D acoustic medium which is highly attenuating material
    spatial_dim = 2
    medium = Acoustic(spatial_dim; ρ=1.0, c=1.0)

    # Plane wave source
    direction = [1.0,1.0]
    source = PlaneSource(medium,direction)

    # a species of circular particles
    ms = MultipleScattering # in case you load a package which has it's own definition of Circle.
    particle_medium = Acoustic(spatial_dim; ρ=0.1, c=0.1)
    specie = Specie(Particle(particle_medium, ms.Circle(0.8)); volume_fraction = 0.25)

    # a halfspace with outward normal [-1.0,0.0] and filled with specie
    material = Material(Halfspace([-1.0,0.0]), specie)

    ω = 1.0
    θin = 0.3 # angle made between incident wave direction and the material boundary

    # calculate the 4 effective waves with the smallest attenuation
    wave_effs = WaveModes(ω, source, material;
          num_wavenumbers = 5,
          θin = θin,
          extinction_rescale = false
    )

    # using Plots
    # plot(wave_effs[1:5], markeralpha = 1.0)
    # savefig("wave_effs.png")
```
![Many effective wavenumbers](wave_effs.png)

## Determine the amplitude/eignvectors of the wavemodes

To determine the amplitude of these effective waves we use the Matching method.
```julia
    match_wave = MatchPlaneWaveMode(ω, source, material;
      θin = θin,
      wave_effs = wave_effs[1:5],
      max_size=80)

    # the error of the match is used as a proxy to determine if the method converged. That is, the match_error is the difference between a sum of the effective waves and a discrete solution.
    match_error(match_wave, material.shape) < 1e-5

    # plot(match_wave, ylim = (-0.01,0.04), xlim = (2.0,10.0),
    #     xlab = "x", hankel_indexes=0:1
    # )
    # savefig("match_wave.png")
```
![The matched wave field](match_wave.png)

We can compare this results against a purely numerical method, which requires a much larger mesh.
```julia
  discrete_wave = DiscretePlaneWaveMode(ω, source, material; θin = θin, max_size=700)

  plot(discrete_wave, xlab = "x", hankel_indexes=0:1)
  plot!(match_wave, seriestype=:line, xlims = (0.0,6.0),
      match_region=false, hankel_indexes=0:1
  )
  savefig("compare_match_wave.png")
```
The two fields perfectly overlap, as the matched method is an exact method. Below the line is from the matched method, while the scatter points are from the purely numerical method.
![Compare matched and discrete field](compare_match_wave.png)

## Calculate the reflection coefficient
```julia

R = reflection_coefficient(ω, match_wave, source, material)
```
