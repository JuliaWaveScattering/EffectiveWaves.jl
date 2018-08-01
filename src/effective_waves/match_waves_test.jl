using EffectiveWaves
using OffsetArrays

using Plots; pyplot()

include("../average_waves/integral_form.jl")

# physical parameters
medium = Medium(1.0,1.0+0.0im)
specie = Specie(ρ=0.5, r=0.5, c=0.2, volfrac=0.1)

k = 0.8;
θin = 0.0
ω = real(k*medium.c)
ho = 1

T = Float64
tol = 1e-7
X = 0.0:0.05:8.
XJ = length(X)

XL = XJ - 8
# XL = Int(round(length(X)/2))
X_match = X[XL]

# Calculate discretised average wave directly from governing equations
avg_wave = AverageWave(ω, medium, specie; hankel_order=ho, X=X)

# From effective wave theory
k_effs = wavenumbers(ω, medium, [specie];
    mesh_size = 0.2, mesh_points = 20,
    tol = tol/10, hankel_order = ho
)
scatter(real.(k_effs),imag.(k_effs), title = "Effective wavenumbers", xlab = "Re k_eff", ylab = "Im k_eff")

wave_effs = [
    EffectiveWave(ω, k_eff, medium, [specie];
        hankel_order = ho, extinction_rescale = true, tol=tol*10
    )
for k_eff in k_effs]

# rescale to lessen roundoff errors
for i in eachindex(k_effs)[2:end]
    wave_effs[i].amplitudes = wave_effs[i].amplitudes * exp( -im*k_effs[i]*X_match) / norm(wave_effs[i].amplitudes)
end

species = [specie]
include("match_waves.jl")

(w_vec, q_arr) = extinc_arrays(ω, wave_effs, X_match, X, medium, species; θin = θin)
(L_mat, E_mat, b_eff) = match_arrays(ω, wave_effs, X_match, X, medium, species; θin = θin)

LA = L_mat*avg_wave.amplitudes[:]
qA = sum(q_arr[ind] * avg_wave.amplitudes[ind] for ind in eachindex(avg_wave.amplitudes))

αs = LA + b_eff

# check extinction
    norm((sum(w_vec .* αs) - (qA + im*k^2))/sum(w_vec .* αs)) < 1e-12

# check matching
    # Calculate the discretised wave from these effective wave
    avg_wave_effs = [AverageWave(real(k), wave, X) for wave in wave_effs]
    amplitudes_eff = sum(αs[i] * avg_wave_effs[i].amplitudes[XL:end,:,:] for i in eachindex(avg_wave_effs))
    avg_wave_eff = AverageWave(ho, X[XL:end], amplitudes_eff)

    plot(xlab = "X", ylab = "Re amplitude")
    for n = -ho:ho
        scatter!(avg_wave_eff.X, real.(avg_wave_eff.amplitudes[:,n+ho+1]), lab = "Wave eff hankel $n")
        plot!(X[XL:end], real.(avg_wave.amplitudes[XL:end,n+ho+1]), lab = "Avg wave hankel $n" )
    end
    gui()

    plot(xlab = "X", ylab = "Abs amplitude")
    for n = -ho:ho
        scatter!(avg_wave_eff.X, abs.(avg_wave_eff.amplitudes[:,n+ho+1]), lab = "Wave eff hankel $n")
        plot!(X[XL:end], abs.(avg_wave.amplitudes[XL:end,n+ho+1]), lab = "Avg wave hankel $n" )
    end
    gui()
