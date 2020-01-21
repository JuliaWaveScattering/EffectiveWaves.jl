using EffectiveWaves
using OffsetArrays

using Plots; pyplot()

include("../discrete_wave/integral_form.jl")

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
avg_wave = DiscretePlaneWaveMode(ω, medium, specie; basis_order=ho, X=X)

# From effective wave theory
k_effs = wavenumbers(ω, medium, [specie];
    mesh_size = 0.2, mesh_points = 20,
    tol = tol/10, basis_order = ho
)
scatter(real.(k_effs),imag.(k_effs), title = "Effective wavenumbers", xlab = "Re k_eff", ylab = "Im k_eff")

wave_effs = [
    EffectivePlaneWaveMode(ω, k_eff, medium, [specie];
        basis_order = ho, extinction_rescale = true, tol=tol*10
    )
for k_eff in k_effs]

# rescale to lessen roundoff errors
for i in eachindex(k_effs)[2:end]
    wave_effs[i].amplitudes = wave_effs[i].amplitudes * exp( -im*k_effs[i]*X_match) / norm(wave_effs[i].amplitudes)
end

# Calculate the discretised wave from these effective wave
avg_wave_effs = [DiscretePlaneWaveMode(X, wave_effs[1])]
avg_wave_effs = [[DiscretePlaneWaveMode(X, wave) for wave in wave_effs[2:end]]; avg_wave_effs]

species = [specie]
include("match_waves.jl")

vs = [
    [w.amplitudes[j,n+ho+1,1] for w in avg_wave_effs]
for j = XL:XJ, n = -ho:ho]

invV = inv(sum(vs[inds] * transpose(vs[inds])  for inds in eachindex(vs)))
inv_w = one(T)/(transpose(w_vec)*invV*w_vec)

invVY = hcat(
    [
        invV * [(j < XL) ? zero(Complex{T}) : w.amplitudes[j,n+ho+1,1] for w in avg_wave_effs]
    for j = 1:XJ, n = -ho:ho]...
)

# Yns = [
#     [zeros(Complex{T},length(wave_effs),i_match-1) hcat(vs[:,n+ho+1]...)]
# for n = -ho:ho]
# Y = hcat(Yns...)
# L = invV*Y + (w_vec.*inv_w) * (transpose(q_arr[:]) - transpose(w_vec)*invV*Y)

invVY_mat = [
    invV * [w.amplitudes[j,n+ho+1,1] for w in avg_wave_effs]
for j = 1:XJ, n = -ho:ho];

invVYA = sum(
    invVY_mat[j,n+ho+1] * avg_wave.amplitudes[j,n+ho+1]
for j = XL:XJ, n = -ho:ho)

norm(invVY  * avg_wave.amplitudes[:] - invVYA)/norm(invVYA) < 1e-12

qA = sum(q_arr[ind] * avg_wave.amplitudes[ind] for ind in eachindex(avg_wave.amplitudes))
LA = invVYA + invV * (w_vec.*inv_w) * (qA - transpose(w_vec)*invVYA)

L = invVY + invV * (w_vec.*inv_w) * (transpose(q_arr[:]) - transpose(w_vec)*invVY)
norm(LA - L*avg_wave.amplitudes[:])/(norm(LA)) < 1e-12

αs = LA + (im*k^2*inv_w).*invV*w_vec

# check extinction
norm((sum(w_vec .* αs) - (qA + im*k^2))/sum(w_vec .* αs))

# check matching
amplitudes_eff = sum(αs[i] * avg_wave_effs[i].amplitudes[XL:end,:,:] for i in eachindex(avg_wave_effs))
avg_wave_eff = DiscretePlaneWaveMode(ho, X[XL:end], amplitudes_eff)

plot(xlab = "X", ylab = "Re amplitude")
for n = -ho:ho
    scatter!(avg_wave_eff.x, real.(avg_wave_eff.amplitudes[:,n+ho+1]), lab = "Wave eff hankel $n")
    plot!(X[XL:end], real.(avg_wave.amplitudes[XL:end,n+ho+1]), lab = "Avg wave hankel $n" )
end
gui()

plot(xlab = "X", ylab = "Abs amplitude")
for n = -ho:ho
    scatter!(avg_wave_eff.x, abs.(avg_wave_eff.amplitudes[:,n+ho+1]), lab = "Wave eff hankel $n")
    plot!(X[XL:end], abs.(avg_wave.amplitudes[XL:end,n+ho+1]), lab = "Avg wave hankel $n" )
end
gui()
