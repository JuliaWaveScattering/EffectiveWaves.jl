using EffectiveWaves
using OffsetArrays

include("match_waves.jl")
include("../average_waves/integral_form.jl")

# physical parameters
medium = Medium(1.0,1.0+0.0im)
specie = Specie(ρ=0.5, r=0.5, c=0.2, volfrac=0.1)

k = 0.8;
θin = 0.0
ω = real(k*medium.c)
ho = 1

tol = 1e-5
X = 0.0:0.01:3.

# From effective wave theory
k_effs = wavenumbers(ω, medium, [specie]; mesh_points = 10, tol = tol, hankel_order = ho)
wave_effs = [
    EffectiveWave(ω, k_eff, medium, [specie];
        hankel_order = ho, tol=tol*10, extinction_rescale = false
    )
for k_eff in k_effs[2:end]]
wave_effs = [EffectiveWave(ω, k_effs[1], medium, [specie]; hankel_order = ho, tol=tol*10); wave_effs]

# Calculate the discretised wave from these effective wave
avg_wave_effs = [AverageWave(wave, X; ) for wave in wave_effs]

# Calculate discretised average wave directly from governing equations
avg_wave = AverageWave(ω, medium, specie; hankel_order=ho, X=X)

using Plots; pyplot()

plot(X, abs.(avg_wave.amplitudes[:,2]), xlab = "depth X", lab = "avg. wave")

map(eachindex(k_effs)) do i
    plot!(X, abs.(avg_wave_effs[i].amplitudes[:,2]), xlab = "depth X", lab = "eff. wave k_eff = $(k_effs[i])")
end
gui()

X_match = 0.;
ind_match = findmin(abs.(X .- X_match))[2]
extinc = overlap_arrays(ω, [wave_effs[1]], X_match, X, medium, [specie]; θin = θin)

q = extinc[2][:]
w = extinc[1][1]
a = im*k^2 + dot(q,avg_wave.amplitudes[1:ind_match,:,:][:])
# wave_effs[1].amplitudes = (a/w)*wave_effs[1].amplitudes

αs = map(0.:0.2:2.) do x
    ind_match = findmin(abs.(X .- x))[2]
    extinc = overlap_arrays(ω, [wave_effs[1]], x, X, medium, [specie]; θin = θin)

    q = extinc[2][:]
    w = extinc[1][1]
    a = im*k^2 + dot(q,avg_wave.amplitudes[1:ind_match,:,:][:])
    a/w
end
abs.(αs)



# Calculate the discretised wave from these effective wave
avg_wave_effs = [AverageWave(wave, X, X_match) for wave in wave_effs]

plot(X, abs.(avg_wave.amplitudes[:,2]), xlab = "depth X", lab = "avg. wave")

map(eachindex(k_effs[1:1])) do i
    plot!(avg_wave_effs[i].X, abs.(avg_wave_effs[i].amplitudes[:,2]), xlab = "depth X", lab = "eff. wave k_eff = $(k_effs[i])")
end
gui()
