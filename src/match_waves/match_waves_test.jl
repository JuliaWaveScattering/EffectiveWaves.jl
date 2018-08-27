using EffectiveWaves
using OffsetArrays
using Plots; pyplot()
Plots.scalefontsizes(1.5)

# physical parameters
medium = Medium(1.0,1.0+0.0im)
specie = Specie(ρ=0.5, r=0.4, c=0.2, volfrac=0.1)
specie = Specie(ρ=0.2, r=0.4, c=0.2, volfrac=0.2)

k = 1.;
θin = 0.0
ω = real(k*medium.c)
radius_multiplier = 1.005

kws = ()
T = Float64
tol = 1e-7
# ho = maximum_hankel_order(ω, medium, [specie]; tol=tol)
ho = 2

# for X = 0.0:0.05:12. the relative convergence error for
# Specie(ρ=0.5, r=0.5, c=0.2, volfrac=0.1), k = 0.8, and θin = 0.0,
#  inbetween 6 and 8 is 0.05%
ka12 = k*radius_multiplier*2.0*specie.r
X = 0.0:(ka12/31):15.
# X = 0.0:0.04:12.

# Calculate discretised average wave directly from governing equations
avg_wave_X = AverageWave(ω, medium, specie; hankel_order=ho, θin=θin, X = X,
    radius_multiplier = radius_multiplier)
# x = X./k

# plot(avg_wave_X, hankel_indexes = 0:0, seriestype=:scatter)
# plot!(avg_wave_X, hankel_indexes = 0:0, apply=imag, seriestype=:scatter)

# wave_effs = effective_waves(ω, medium, [specie];
#     hankel_order=ho, tol = tol,  θin = θin,
#     radius_multiplier = radius_multiplier, mesh_points = 6
#     , extinction_rescale = false
#     )
#
# k_effs = [w.k_eff for w in wave_effs]
# scatter(real.(k_effs),imag.(k_effs),
#     title = "Effective wavenumbers",
#     xlab = "Re k_eff", ylab = "Im k_eff",
#     ylims = (0,Inf))
include("match_waves.jl")

match_w = match_effective_waves(ω, medium, specie;
    radius_multiplier = radius_multiplier,
    hankel_order=ho,
    tol = tol, θin = θin)
    # ,wave_effs = wave_effs);

L = match_w.match_index
x = match_w.average_wave.X[L]:0.01:avg_wave_X.X[end]

plot(avg_wave_X, hankel_indexes = 1:2, seriestype=:scatter)
plot!(avg_wave_X, hankel_indexes = 1:2, apply=abs, seriestype=:scatter)

plot!(match_w.average_wave, hankel_indexes = 1:2)
plot!(match_w.average_wave, hankel_indexes = 1:2, apply=abs)

plot!(x, match_w.effective_waves, hankel_indexes = 1:2)
plot!(x, match_w.effective_waves, hankel_indexes = 1:2, apply=abs)

# X = avg_wave_X.X
X0 = 2. # match from X[L]
L = findmin(abs.(X .- X0))[2]
X_max = X[Int(round(0.65*length(X)))] # maximum X considered to have the correct numerical solution
i_max = findmin(abs.(X .- X_max))[2]


# # If we assume the average wave is a sum of plane waves everywhere, then we can estimate which of these effectives waves will have decayed too much already, and therefore discard them
avg_wave_effs = [AverageWave(real(k), wave, X[L:L+1]) for wave in wave_effs]
ref_norm = norm(avg_wave_effs[1].amplitudes[end,:,:])/norm(wave_effs[1].amplitudes)
eff_inds = find(
    norm(avg_wave_effs[i].amplitudes[end,:,:])/norm(wave_effs[i].amplitudes) > ref_norm*tol
for i in eachindex(wave_effs))

k_effs = k_effs[eff_inds]
wave_effs = wave_effs[eff_inds];

J = L + 3*length(k_effs) # double number of points than waves to avoid overfitting
X_match = X[L:J]
avg_wave = AverageWave(ho, X[1:J], avg_wave_X.amplitudes[1:J,:,:])

# rescale to lessen roundoff errors
for i in eachindex(k_effs)[2:end]
    wave_effs[i].amplitudes = wave_effs[i].amplitudes * exp(-im*k_effs[i]*cos(wave_effs[i].θ_eff)*X_match[1]/k) / norm(wave_effs[i].amplitudes)
end

include("match_waves.jl")

Y = match_only_arrays(ω, wave_effs, L, X[1:J], medium, [specie]; θin = θin)
(w_vec, G_arr) = extinc_arrays(ω, wave_effs, L, X[1:J], medium, [specie]; θin = θin)
(LT_mat, E_mat, b_eff) = match_arrays(ω, wave_effs, L, X[1:J], medium, [specie]; θin = θin)

LA = LT_mat*avg_wave.amplitudes[:]
qA = sum(G_arr[ind] * avg_wave.amplitudes[ind] for ind in eachindex(avg_wave.amplitudes))

αs = LA + b_eff
scatter([abs.(exp.(im*k_effs./k)), sign.(real.(k_effs)).*abs.(αs)], lab = ["exp(im*k_eff)" "sign(Re k_eff) * abs αs"])


# check extinction
    norm((sum(w_vec .* αs) - (qA + im*k^2))/sum(w_vec .* αs)) < 1e-10

# For only one effective wave
    (w_vec, q_arr) = extinc_arrays(ω, wave_effs[1:1], L, X[1:J], medium, [specie]; θin = θin);
    qA = sum(q_arr[ind] * avg_wave.amplitudes[ind] for ind in eachindex(avg_wave.amplitudes));

    α1 = (qA + im*k^2)/w_vec[1]

# check matching
    αs_match = Y*avg_wave.amplitudes[:]
    # Calculate the discretised wave from these effective wave
    avg_wave_effs = [AverageWave(real(k), wave, X[L:i_max]) for wave in wave_effs]

    amplitudes_eff = sum(αs[i] * avg_wave_effs[i].amplitudes[:,:,:] for i in eachindex(avg_wave_effs))
    avg_wave_eff = AverageWave(ho, avg_wave_effs[1].X, amplitudes_eff)

    amplitudes_eff = sum(αs_match[i] * avg_wave_effs[i].amplitudes[:,:,:] for i in eachindex(avg_wave_effs))
    avg_wave_eff_match = AverageWave(ho, avg_wave_effs[1].X, amplitudes_eff)

    ps = map(0:ho) do n
        max_w = maximum(real.(avg_wave.amplitudes[L:end,n+ho+1]))
        min_w = minimum(real.(avg_wave.amplitudes[L:end,n+ho+1]))
        max_w = (max_w > 0) ? 1.2*max_w : 0.0
        min_w = (min_w < 0) ? 1.2*min_w : 0.0
        plot(x -> x,y->max_w,X[L],X[J],line=0,fill=(0,:orange), fillalpha=0.4, lab="match region")
        plot!(x -> x,y->min_w,X[L],X[J],line=0,fill=(0,:orange), fillalpha=0.4, lab="")
        plot!(xlab = "X", ylab = "Re amplitude", title="Hankel order $n")
        plot!(avg_wave_eff.X, real.(avg_wave_eff.amplitudes[:,n+ho+1]), lab = "Wave eff", linewidth=2)
        # plot!(avg_wave_eff_match.X, real.(avg_wave_eff_match.amplitudes[:,n+ho+1]), lab = "Wave eff match", linewidth=2)
        scatter!(X[L:i_max], real.(avg_wave_X.amplitudes[L:i_max,n+ho+1]), lab = "Avg wave", linewidth =2 )
        plot!(avg_wave_effs[1].X, real.(α1*avg_wave_effs[1].amplitudes[:,n+ho+1]), lab = "Wave eff 1", linewidth =2, linestyle=:dash)
    end
    plot(ps...)











    plot(xlab = "X", ylab = "Abs amplitude")
    for n = -ho:ho
        scatter!(avg_wave_eff.X, abs.(avg_wave_eff.amplitudes[:,n+ho+1]), lab = "Wave eff hankel $n")
        plot!(X[L:end], abs.(avg_wave_X.amplitudes[L:end,n+ho+1]), lab = "Avg wave hankel $n" )
    end
    gui()
