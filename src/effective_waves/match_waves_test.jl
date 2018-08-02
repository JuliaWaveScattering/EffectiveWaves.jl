using EffectiveWaves
using OffsetArrays

include("../average_waves/integral_form.jl")

# physical parameters
medium = Medium(1.0,1.0+0.0im)
specie = Specie(ρ=0.5, r=0.5, c=0.2, volfrac=0.1)

k = 0.8;
θin = 0.0
ω = real(k*medium.c)

T = Float64
tol = 1e-7
# ho = 3
ho = maximum_hankel_order(ω, medium, [specie]; tol=tol*1e4)

X = 0.0:0.05:12.
# XJ = length(X)
XJ = Int(round(length(X)/2 + length(X)/10))
XL = Int(round(length(X)/2 - length(X)/10))
XJ = Int(round(length(X)/2 + 20))
XL = Int(round(length(X)/2 - 20))
Xend =  Int(round(length(X)/2 + length(X)/4))
X_match = X[XL:XJ]

# Calculate discretised average wave directly from governing equations
avg_wave_X = AverageWave(ω, medium, specie; hankel_order=ho, X=X)
avg_wave = AverageWave(ho, X[1:XJ], avg_wave_X.amplitudes[1:XJ,:,:])

# From effective wave theory
k_effs = wavenumbers(ω, medium, [specie];
    mesh_size = 0.2, mesh_points = 20,
    tol = tol/10, hankel_order = ho
)

using Plots; pyplot()

scatter(real.(k_effs),imag.(k_effs), title = "Effective wavenumbers", xlab = "Re k_eff", ylab = "Im k_eff")

# If we assume the average wave is a sum of plane waves every where, then we can estimate which of these effectives waves will have decayed too much already, and therefore discard them
decay_effs = abs.(exp.((im*X_match[1]/k).*k_effs))
eff_inds = find(decay_effs .> 1e-8)
k_effs = k_effs[eff_inds]

wave_effs = [
    EffectiveWave(ω, k_eff, medium, [specie];
        hankel_order = ho, extinction_rescale = true, tol=tol*10
    )
for k_eff in k_effs]

# rescale to lessen roundoff errors
for i in eachindex(k_effs)[2:end]
    wave_effs[i].amplitudes = wave_effs[i].amplitudes * exp(-im*k_effs[i]*X_match[1]/k) / norm(wave_effs[i].amplitudes)
end

include("match_waves.jl")

(w_vec, q_arr) = extinc_arrays(ω, wave_effs, XL, X[1:XJ], medium, [specie]; θin = θin)
(L_mat, E_mat, b_eff) = match_arrays(ω, wave_effs, XL, X[1:XJ], medium, [specie]; θin = θin)

LA = L_mat*avg_wave.amplitudes[:]
qA = sum(q_arr[ind] * avg_wave.amplitudes[ind] for ind in eachindex(avg_wave.amplitudes))

αs = LA + b_eff

# check extinction
    norm((sum(w_vec .* αs) - (qA + im*k^2))/sum(w_vec .* αs)) < 1e-12

# check matching
    # Calculate the discretised wave from these effective wave
    avg_wave_effs = [AverageWave(real(k), wave, X[XL:Xend]) for wave in wave_effs]
    amplitudes_eff = sum(αs[i] * avg_wave_effs[i].amplitudes[:,:,:] for i in eachindex(avg_wave_effs))
    avg_wave_eff = AverageWave(ho, avg_wave_effs[1].X, amplitudes_eff)

    ps = map(0:ho) do n
        max_w = maximum(real.(avg_wave.amplitudes[XL:end,n+ho+1]))
        min_w = minimum(real.(avg_wave.amplitudes[XL:end,n+ho+1]))
        max_w = (max_w > 0) ? 1.2*max_w : 0.8*max_w
        min_w = (min_w < 0) ? 1.2*min_w : 0.8*min_w
        plot(x -> x,y->max_w,X[XL],X[XJ],line=0,fill=(0,:orange), fillalpha=0.4, lab="match region")
        plot!(x -> x,y->min_w,X[XL],X[XJ],line=0,fill=(0,:orange), fillalpha=0.4, lab="")
        plot!(xlab = "X", ylab = "Re amplitude", title="Hankel order $n")
        plot!(avg_wave_eff.X, real.(avg_wave_eff.amplitudes[:,n+ho+1]), lab = "Wave eff")
        scatter!(X[XL:Xend], real.(avg_wave_X.amplitudes[XL:Xend,n+ho+1]), lab = "Avg wave", linewidth =2 )
        plot!(avg_wave_effs[1].X, real.(avg_wave_effs[1].amplitudes[:,n+ho+1]), lab = "Wave eff 1", linewidth =2, linestyle=:dash)
    end
    plot(ps...)











    plot(xlab = "X", ylab = "Abs amplitude")
    for n = -ho:ho
        scatter!(avg_wave_eff.X, abs.(avg_wave_eff.amplitudes[:,n+ho+1]), lab = "Wave eff hankel $n")
        plot!(X[XL:end], abs.(avg_wave_X.amplitudes[XL:end,n+ho+1]), lab = "Avg wave hankel $n" )
    end
    gui()
