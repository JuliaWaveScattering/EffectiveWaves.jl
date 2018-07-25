# physical parameters
medium = Medium(1.0,1.0+0.0im)
specie = Specie(ρ=0.5, r=0.5, c=0.2, volfrac=0.1)

k = 0.8;
ω = real(k*medium.c)
ho = 1

tol = 1e-5
x = 0.0:0.01:2.

# From effective wave theory
k_effs = wavenumbers(ω, medium, [specie]; mesh_points = 10, tol = tol, hankel_order = ho)
wave_effs = [
    EffectiveWave(ω, k_eff, medium, [specie];
        hankel_order = ho, tol=tol*10, extinction_rescale = false
    )
for k_eff in k_effs[2:end]]
wave_effs = [EffectiveWave(ω, k_effs[1], medium, [specie]; hankel_order = ho, tol=tol*10); wave_effs]

# Calculate the discretised wave from these effective wave
avg_wave_effs = [AverageWave(ω, x, wave) for wave in wave_effs]

# Calculate discretised average wave directly from governing equations
avg_wave = AverageWave(ω, medium, specie; hankel_order=ho, x=x)

using Plots; pyplot()

plot(x, abs.(avg_wave.amplitudes[:,2]), xlab = "depth x", lab = "avg. wave")

map(eachindex(k_effs)) do i
    plot!(x, abs.(avg_wave_effs[i].amplitudes[:,2]), xlab = "depth x", lab = "eff. wave k_eff = $(k_effs[i])")
end
gui()
