using EffectiveWaves
using JLD2
using FileIO
# using OffsetArrays
# using Plots; pyplot()
# Plots.scalefontsizes(1.5)

# physical parameters
medium = Medium(1.0,1.0+0.0im)
specie = Specie(ρ=0.1, r=0.4, c=0.1, volfrac=0.04)
specie = Specie(ρ=0.2, r=0.4, c=0.2, volfrac=0.2)
# specie = Specie(ρ=20.0, r=0.4, c=20.0, volfrac=0.10)
species = [specie]

k = 1.0
# k = 0.4;
θin = 0.4
# θin = 0.0
ω = real(k*medium.c)
radius_multiplier = 1.005

kws = ()
T = Float64
tol = 1e-7
ho = maximum_basis_order(ω, medium, [specie]; tol=tol) -1
basis_order = ho

# for X = 0.0:0.05:12. the relative convergence error for
# Specie(ρ=0.5, r=0.5, c=0.2, volfrac=0.1), k = 0.8, and θin = 0.0,
#  inbetween 6 and 8 is 0.05%
a12 = radius_multiplier*2.0*specie.r
num_wavenumbers = 90;
num_wavenumbers = 40;
# X = 0.0:(a12k/31):15.
# X = 0.0:0.04:12.


wave_effs = effective_waves(ω, medium, [specie];
    basis_order=ho, tol = 1e-7,  θin = θin, verbose = true,
    num_wavenumbers = num_wavenumbers,
    apply_meshing = true,
    radius_multiplier = radius_multiplier,
    mesh_refine = 0.18, mesh_points = 10, mesh_size=2.0,
    extinction_rescale = false
)

# k_effs = [w.k_eff for w in wave_effs];
# mesh_ks = wavenumbers_mesh(ω, k_effs, medium, [specie];
#     basis_order=ho, tol = 1e-7,  θin = θin, verbose = true,
#     num_wavenumbers = num_wavenumbers,
#     radius_multiplier = radius_multiplier, mesh_refine = 0.2
#     );

using Plots; pyplot()
plot(wave_effs, color=:black)
plot!(wave_effs[1:num_wavenumbers], markercolor=:white)

collect(x_mesh_match(wave_effs[1:num_wavenumbers])[2])

# Calculate discretised average wave directly from governing equations
L, x = x_mesh_match(wave_effs[1:num_wavenumbers]; a12 = a12, max_size = 220, tol=1e-9);

x = x_mesh(wave_effs[1], wave_effs[2]; a12 = a12, max_size = 1600, tol=1e-8)

avg_wave_x = AverageWave(ω, medium, specie;
    basis_order=ho, θin=θin, x = x,
    wave_effs = wave_effs[1:num_wavenumbers],
    radius_multiplier = radius_multiplier, tol = 1e-8)
X = avg_wave_x.x.*k
# x = X./k

scatter(avg_wave_x)
# scatter!(avg_wave_x, hankel_indexes = 0:ho, apply=abs)

L, x = x_mesh_match(wave_effs[1:num_wavenumbers]; a12 = a12, max_size = 220, tol=1e-9);
x_short = x[1:Int(round(length(x)/3))]

match_w = MatchWave(ω, medium, specie;
    x = x_short, L_match = 20,
    radius_multiplier = radius_multiplier,
    basis_order=ho,
    tol = 1e-7, θin = θin, wave_effs = wave_effs[1:num_wavenumbers]);

    ω = real(k*medium.c)
    radius_multiplier = 1.005

    kws = ()
    T = Float64
    tol = 1e-7
    ho = maximum_basis_order(ω, medium, [specie]; tol=tol) -1
    basis_order = ho

plot(match_w, hankel_indexes = 0:1)
plot(avg_wave_x, hankel_indexes = 0:1, seriestype=:scatter, xlims = (0.,5.))
plot!(match_w, hankel_indexes = 0:1, seriestype=:line, xlims = (0.,5.))

avg_wave_small = AverageWave(ω, medium, specie; basis_order=ho, θin=θin, x=match_w.discrete_wave.x,
    radius_multiplier = radius_multiplier, tol = 1e-6, max_size=500)
plot!(avg_wave_small, seriestype=:line, linestyle=:dash)

# X = avg_wave_x.x
X0 = 2. # match from X[L]
L = findmin(abs.(X .- X0))[2]
X_max = X[Int(round(0.65*length(X)))] # maximum X considered to have the correct numerical solution
i_max = findmin(abs.(X .- X_max))[2]


# # If we assume the average wave is a sum of plane waves everywhere, then we can estimate which of these effectives waves will have decayed too much already, and therefore discard them
avg_wave_effs = [AverageWave(X[L:L+1], wave) for wave in wave_effs]
ref_norm = norm(avg_wave_effs[1].amplitudes[end,:,:])/norm(wave_effs[1].amplitudes)
eff_inds = find(
    norm(avg_wave_effs[i].amplitudes[end,:,:])/norm(wave_effs[i].amplitudes) > ref_norm*tol
for i in eachindex(wave_effs))

k_effs = k_effs[eff_inds]
wave_effs = wave_effs[eff_inds];

J = L + 3*length(k_effs) # double number of points than waves to avoid overfitting
X_match = X[L:J]
avg_wave = AverageWave(ho, X[1:J], avg_wave_x.amplitudes[1:J,:,:])

# rescale to lessen roundoff errors
for i in eachindex(k_effs)[2:end]
    wave_effs[i].amplitudes = wave_effs[i].amplitudes * exp(-im*k_effs[i]*cos(wave_effs[i].θ_eff)*X_match[1]/k) / norm(wave_effs[i].amplitudes)
end

include("match_and_extinction.jl")

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
    avg_wave_effs = [AverageWave(X[L:i_max], wave) for wave in wave_effs]

    amplitudes_eff = sum(αs[i] * avg_wave_effs[i].amplitudes[:,:,:] for i in eachindex(avg_wave_effs))
    avg_wave_eff = AverageWave(ho, avg_wave_effs[1].x, amplitudes_eff)

    amplitudes_eff = sum(αs_match[i] * avg_wave_effs[i].amplitudes[:,:,:] for i in eachindex(avg_wave_effs))
    avg_wave_eff_match = AverageWave(ho, avg_wave_effs[1].x, amplitudes_eff)

    ps = map(0:ho) do n
        max_w = maximum(real.(avg_wave.amplitudes[L:end,n+ho+1]))
        min_w = minimum(real.(avg_wave.amplitudes[L:end,n+ho+1]))
        max_w = (max_w > 0) ? 1.2*max_w : 0.0
        min_w = (min_w < 0) ? 1.2*min_w : 0.0
        plot(x -> x,y->max_w,X[L],X[J],line=0,fill=(0,:orange), fillalpha=0.4, lab="match region")
        plot!(x -> x,y->min_w,X[L],X[J],line=0,fill=(0,:orange), fillalpha=0.4, lab="")
        plot!(xlab = "X", ylab = "Re amplitude", title="Hankel order $n")
        plot!(avg_wave_eff.x, real.(avg_wave_eff.amplitudes[:,n+ho+1]), lab = "Wave eff", linewidth=2)
        # plot!(avg_wave_eff_match.x, real.(avg_wave_eff_match.amplitudes[:,n+ho+1]), lab = "Wave eff match", linewidth=2)
        scatter!(X[L:i_max], real.(avg_wave_x.amplitudes[L:i_max,n+ho+1]), lab = "Avg wave", linewidth =2 )
        plot!(avg_wave_effs[1].x, real.(α1*avg_wave_effs[1].amplitudes[:,n+ho+1]), lab = "Wave eff 1", linewidth =2, linestyle=:dash)
    end
    plot(ps...)




    plot(xlab = "X", ylab = "Abs amplitude")
    for n = -ho:ho
        scatter!(avg_wave_eff.x, abs.(avg_wave_eff.amplitudes[:,n+ho+1]), lab = "Wave eff hankel $n")
        plot!(X[L:end], abs.(avg_wave_x.amplitudes[L:end,n+ho+1]), lab = "Avg wave hankel $n" )
    end
    gui()
