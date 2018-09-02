using EffectiveWaves

medium = Medium(1.0,1.0+0.0im)
# Large weak scatterers with low volume fraciton
specie = Specie(ρ=0.5,r=0.5, c=0.2, volfrac=0.1)
ω = 0.5
k = real(ω/medium.c)
hankel_order = 1
tol = 1e-6

k_effs = wavenumbers(ω, medium, [specie]; tol = tol, mesh_points = 10, hankel_order=hankel_order)
wave_effs = [EffectiveWave(ω, k_eff, medium, [specie]; tol = tol/10., extinction_rescale=false) for k_eff in k_effs]

# generate an averaged wave over x from these effective wave
x = 0.:0.01:1.
avg_wave_effs = [AverageWave(x, wave) for wave in wave_effs]

# calculate the numerical averaged wave over x
avg_wave = AverageWave(ω, medium, specie; x = x, hankel_order=hankel_order)

using Plots; pyplot()

# plot(x,[real.(avg_wave.amplitudes[:,1]),real.(avg_wave.amplitudes[:,2])], lab =["avg. hankel -1","avg. hankel 0"])
f = identity
plot(x, f.(abs.(avg_wave.amplitudes[:,2])), lab = "avg. ho =0")

map(eachindex(k_effs)) do i
    plot!(x, f.(abs.(avg_wave_effs[i].amplitudes[:,2])), lab = "k_eff = $(k_effs[i]), ho = 0")
end
gui()
