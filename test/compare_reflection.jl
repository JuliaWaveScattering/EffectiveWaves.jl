using EffectiveWaves


Rs = reflection_coefficient(ωs, k_effs, medium, species)
Rs2 = reflection_coefficient(ωs, medium, species)

@test norm(Rs - Rs2) < 1e-8*norm(Rs)

R_low = reflection_coefficient_halfspace(medium, eff_medium)
Rs_low = reflection_coefficient(ωs, k_eff_lows, medium, species)
Rs = reflection_coefficient(ωs, k_effs, medium, species)

norm(Rs_low[1] - R_low)/norm(R_low)


ωs = collect(0.01:0.01:0.05)
# strong scatterers
species = [
    Specie(ρ=0.1,r=0.1, c=0.02, volfrac=0.05),
    Specie(ρ=3., r=0.2, c=0.01, volfrac=0.04)
]

eff_medium = effective_medium(medium, species)
k_eff_lows = ωs./eff_medium.c

k_eff_φs = wavenumber_low_volfrac(ωs, medium, species)
k_effs = wavenumber(ωs, medium, species)

R_low = reflection_coefficient_halfspace(medium, eff_medium)
Rs_low = reflection_coefficient(ωs, k_eff_lows, medium, species)
Rs = reflection_coefficient(ωs, k_effs, medium, species)

norm(Rs_low[1] - R_low)/norm(R_low)
