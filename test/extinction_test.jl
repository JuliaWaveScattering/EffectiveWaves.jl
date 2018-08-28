# test extinction of the incident wave

using EffectiveWaves
using OffsetArrays

# using Plots; pyplot()

include("../src/average_waves/integral_form.jl")
include("../src/effective_waves/match_waves.jl")


# physical parameters
medium = Medium(1.0,1.0+0.0im)
specie = Specie(ρ=0.5, r=0.5, c=0.2, volfrac=0.1)

k = 0.8;
θin = 0.0
ω = real(k*medium.c)
ho = 1

tol = 1e-7
X = 0.0:0.05:8.


function qA_eff(w::EffectiveWave, Xi)
    (2*k*specie.num_density/cos(θin)) *
    (exp(im*Xi*(w.k_eff*cos(w.θ_eff)/k - cos(θin))) - 1.)/(w.k_eff*cos(w.θ_eff) - k*cos(θin)) *
    sum(
        Zn(ω,specie,medium,n)*w.amplitudes[n+ho+1]*exp(im*n*(θin - w.θ_eff))
    for n = -ho:ho)
end

# Calculate discretised average wave directly from governing equations
avg_wave = AverageWave(ω, medium, specie; hankel_order=ho, X=X)

# From effective wave theory
k_effs = wavenumbers(ω, medium, [specie]; mesh_points = 10, tol = tol, hankel_order = ho)
wave_effs = [
    EffectiveWave(ω, k_eff, medium, [specie];
        hankel_order = ho, tol=tol*10, extinction_rescale = false
    )
for k_eff in k_effs[2:end]]
wave_effs = [
    EffectiveWave(ω, k_effs[1], medium, [specie];
        hankel_order = ho, extinction_rescale = true, tol=tol*10)
; wave_effs]

# Calculate the discretised wave from these effective wave
avg_wave_effs = [AverageWave(k, wave, X) for wave in wave_effs]


# the as should tend to
inds = Int.(round.(collect(linspace(1,length(X),10))))
das = map(inds) do i
    extinc = extinc_arrays(ω, [wave_effs[1]], X[i], X, medium, [specie]; θin = θin)

    q = extinc[2]
    w = extinc[1][1]
    qA_eff1 = sum(q[ind]*avg_wave_effs[1].amplitudes[ind] for ind in eachindex(q))
    qA = sum(q[ind]*avg_wave.amplitudes[ind] for ind in eachindex(q))

    a = im*k^2 + qA
    # a = im*k^2 +  qA_eff1
    [im*k^2 + qA,  im*k^2 + qA_eff1, im*k^2 + qA_eff(wave_effs[1],X[i]), (a/w)]
end
dqA = [abs(a[1]) for a in das];
dqA_eff1 = [abs(a[2]) for a in das];
dqA_eff = [abs(a[3]) for a in das];
αs = [a[4] for a in das];
# plot(X[inds], [dqA, dqA_eff1, dqA_eff, abs.(αs), abs.(exp.(im*wave_effs[1].k_eff.*X[inds]))],
#     lab = ["extinc Inc" "extinc Inc eff Num" "extinc Inc eff" "aw" "wave decay"], xlab = "X_bar"
# )

norm(dqA - dqA_eff1)/norm(dqA) < 0.05

# Show that the numerical solution becomes more like the least attenuating wave away from both intefaces
    i = Int(round(length(inds)/2))+1
    norm(avg_wave.amplitudes[inds[i]:end,:,1] - αs[i]*avg_wave_effs[1].amplitudes[inds[i]:end,:,1]) /
    norm(avg_wave.amplitudes[inds[i]:end,:,1]) < 0.3
    # wave_diffs = [
    #     norm(avg_wave.amplitudes[inds[i]:end,:,1] - αs[i]*avg_wave_effs[1].amplitudes[inds[i]:end,:,1]) /
    #     norm(avg_wave.amplitudes[inds[i]:end,:,1])
    # for i in eachindex(αs)]
    # plot(X[inds], wave_diffs,
    #     xlab = "depth X", lab = "eff. wave[1] diff")

# n = 1
# plot(avg_wave.x, 0.0.*avg_wave.x, xlab = "depth X", lab = "abs avg wave ")
# map(eachindex(αs)) do i
#     plot!(X[inds[i]:end], abs.(avg_wave.amplitudes[inds[i]:end,n+1+ho] - αs[i]*avg_wave_effs[1].amplitudes[inds[i]:end,n+1+ho]),
#         xlab = "depth X", lab = "eff. wave α = $(round(100*αs[i])/100)")
# end
# gui()
#
# # Compare all effective waves
# plot(X, abs.(avg_wave.amplitudes[:,n+1+ho]), xlab = "depth X", lab = "avg. wave")
#
# map(eachindex(k_effs)) do i
#     plot!(avg_wave_effs[i].x, abs.(avg_wave_effs[i].amplitudes[:,n+1+ho]), xlab = "depth X", lab = "eff. wave k_eff = $(k_effs[i])")
# end
# gui()
