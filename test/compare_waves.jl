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
hankel_order = 3
radius_multiplier = 1.005

incident_medium = medium
θ_inc = 0.2
ω = 0.01
k = ω/incident_medium.c
ωs = collect(linspace(ω,40.0,60))

(β_eff,ρ_eff) = effective_material_properties(incident_medium, species)
k_eff_lows = ωs.*sqrt(ρ_eff/β_eff)

# effective_wavenumber(ω, medium, species)
k_eff_φs = wavenumber_low_volfrac(ωs, medium, species)

# (k_eff,θ_eff,As) = transmitted_planewave(ω, incident_medium, species; θ_inc = θ_inc)
planewaves = [transmitted_planewave(ω, incident_medium, species; θ_inc = θ_inc) for ω in ωs]
k_effs = [w[1] for w in planewaves];
θ_effs = [w[2] for w in planewaves];
As_vec = [w[3] for w in planewaves];

r = mean(s.r for s in species)

plot(r.*ωs./real.(medium.c), [real.(k_effs), imag.(k_effs), k_eff_lows, real(k_eff_φs), imag(k_eff_φs)],
    xlab = "ak",
    labels = ["Re k*" "Im k*" "k low ω" "Re k low φ" "Im k low φ"],
    linewidth=2,
    line =(2,[:solid :solid :dash :dot :dot])
)

## Strong scatterers
species = [
    Specie(ρ=0.8,r=0.1, c=0.2, volfrac=0.1),
    Specie(ρ=0.2, r=0.2, c=0.1, volfrac=0.1)
]

ωs = collect(linspace(ω,40.0,60))

(β_eff,ρ_eff) = effective_material_properties(incident_medium, species)
k_eff_lows = ωs.*sqrt(ρ_eff/β_eff)

# effective_wavenumber(ω, medium, species)
k_eff_φs = wavenumber_low_volfrac(ωs, medium, species)

# (k_eff,θ_eff,As) = transmitted_planewave(ω, incident_medium, species; θ_inc = θ_inc)
planewaves = [transmitted_planewave(ω, incident_medium, species; θ_inc = θ_inc) for ω in ωs]
k_effs = [w[1] for w in planewaves];
θ_effs = [w[2] for w in planewaves];
As_vec = [w[3] for w in planewaves];

r = mean(s.r for s in species)

plot(r.*ωs./real.(medium.c), [real.(k_effs), imag.(k_effs), k_eff_lows, real(k_eff_φs), imag(k_eff_φs)],
    xlab = "ak",
    labels = ["Re k*" "Im k*" "k low ω" "Re k low φ" "Im k low φ"],
    linewidth=2,
    line =(2,[:solid :solid :dash :dot :dot])
)


reflect_medium = Medium(ρ=ρ_eff, c=sqrt(β_eff/ρ_eff))
R_low_ω = reflection_coefficient_halfspace(incident_medium, reflect_medium; θ_inc = θ_inc)
Rs = [reflection_coefficient(ω, incident_medium, species; θ_inc = θ_inc)
