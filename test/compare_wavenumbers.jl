using EffectiveWaves, Memoize, SpecialFunctions, Optim, BlackBoxOptim, NLsolve

Maxtime=100.
T=Float64


"Derivative of Hankel function of the first kind"
function diffhankelh1(n,z)
  if n!=0
    0.5*(hankelh1(-1 + n, z) - hankelh1(1 + n, z))
  else
    - hankelh1(1, z)
  end
end

"Derivative of Bessel function of first kind"
function diffbesselj(n,z)
  if n!=0
    0.5*(besselj(-1 + n, z) - besselj(1 + n, z))
  else
    - besselj(1, z)
  end
end

## Weak scatterers
species = [
    Specie(ρ=10.,r=0.1, c=12., volfrac=0.1),
    Specie(ρ=3., r=0.2, c=2.0, volfrac=0.1)
]
# background medium
medium = Medium(1.0,1.0+0.0im)
hankel_order = :auto
radius_multiplier = 1.005

incident_medium = medium
ω = 0.01
k = ω/incident_medium.c
ωs = collect(linspace(ω,40.0,40))

eff_medium = effective_medium(incident_medium, species)
k_eff_lows = ωs./eff_medium.c

# effective_wavenumber(ω, medium, species)
k_eff_φs = wavenumber_low_volfrac(ωs, medium, species)

# (k_eff,θ_eff,As) = transmitted_planewave(ω, incident_medium, species; θin = θin)
# k_effs = [w[1] for w in planewaves];
k_effs = [wavenumber(ω, incident_medium, species) for ω in ωs];
r = mean(s.r for s in species)

plot(r.*ωs./real.(medium.c), [real.(k_effs), imag.(k_effs), k_eff_lows, real(k_eff_φs), imag(k_eff_φs)],
    xlab = "ak", title = "weak scatterers",
    labels = ["Re k*" "Im k*" "k low ω" "Re k low φ" "Im k low φ"],
    linewidth=2,
    line =(2,[:solid :solid :dash :dot :dot])
)

## Strong scatterers
species = [
    Specie(ρ=0.8, r=0.1, c=0.2, volfrac=0.1),
    Specie(ρ=0.2, r=0.2, c=0.1, volfrac=0.1)
]

eff_medium = effective_medium(incident_medium, species)
k_eff_lows = ωs./eff_medium.c

k_eff_φs = wavenumber_low_volfrac(ωs, medium, species)
k_effs = [wavenumber(ω, incident_medium, species) for ω in ωs]
r = mean(s.r for s in species)

plot(r.*ωs./real.(medium.c), [real.(k_effs), imag.(k_effs), k_eff_lows, real(k_eff_φs), imag(k_eff_φs)],
    xlab = "ak", title = "strong scatterers",
    labels = ["Re k*" "Im k*" "k low ω" "Re k low φ" "Im k low φ"],
    linewidth=2,
    line =(2,[:solid :solid :dash :dot :dot])
)

# This example does not work..

## Super strong scatterers
species = [
    Specie(ρ=0.08, r=0.001, c=0.002, volfrac=0.1),
    Specie(ρ=0.02, r=0.2, c=0.01, volfrac=0.1)
]

k_eff_φs = wavenumber_low_volfrac(ωs, medium, species);
k_effs = [wavenumber(ω, incident_medium, species) for ω in ωs];

r = mean(s.r for s in species)

plot(r.*ωs./real.(medium.c), [real.(k_effs), imag.(k_effs), real(k_eff_φs), imag(k_eff_φs)],
    xlab = "ak", title = "super strong scatterers",
    labels = ["Re k*" "Im k*" "Re k low φ" "Im k low φ"],
    # ylims = (-0.1,k_eff_lows),
    linewidth=2,
    line =(2,[:solid :solid :dash :dot :dot])
)
