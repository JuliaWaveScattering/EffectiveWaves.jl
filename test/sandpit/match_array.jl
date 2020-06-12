import Base.Test: @testset, @test, @test_throws
using OffsetArrays
using EffectiveWaves

include("match_and_extinction.jl")

medium = Medium(1.0,1.0+0.0im)
specie = Specie(ρ=0.2, r=0.4, c=0.2, volfrac=0.2)

k = 1.;
θin = 0.0
ω = real(k*medium.c)
radius_multiplier = 1.005

kws = ()
T = Float64
tol = 1e-5
ho = 2
a12k = k*radius_multiplier*2.0*specie.r

wave_effs = WaveModes(ω, medium, [specie];
    basis_order=ho, tol = 100*tol,  θin = θin,
    radius_multiplier = radius_multiplier, mesh_points = 6
    , extinction_rescale = false
    )
k_effs = [w.k_eff/k for w in wave_effs]

L, X = x_mesh_match(wave_effs; tol = tol/10, a12=a12k)

avg_wave_effs = [DiscretePlaneWaveMode(X[L:L+1], wave) for wave in wave_effs]
for i in eachindex(k_effs)
    wave_effs[i].amplitudes = wave_effs[i].amplitudes / norm(avg_wave_effs[i].amplitudes[1,:,1])
end

Ls = Int.(round.(LinRange(1,length(X) - length(wave_effs)/2,20)))

errs = map(Ls) do l
    # α = YT*A
    YT = match_only_arrays(ω, wave_effs[1:end-1], l, X, medium, [specie]; θin = θin);
    # the most attenuating wave is added as an error!
    amps = sum(DiscretePlaneWaveMode(X, w).amplitudes for w in wave_effs[1:2:(end-1)]);
    # plus a random distributed error
    amps = amps + tol*rand(size(amps)...)/mean(abs.(amps))

    α = zeros(length(wave_effs)-1)
    α[1:2:end] .= 1.
    YT*amps[:]
    mean(abs.(YT*amps[:] - α))
end
@test maximum(errs) < 10*tol
@test maximum(errs[1:10]) < tol/100
@test issorted(errs[5:end])

# Below is to test the approximation to the E integral. But now match_arrays.jl returns E+R.., so code has just been pasted below
J, X = x_mesh_match(wave_effs; tol = tol/200, a12=a12k)
XJ = 0.:(X[2] - X[1]):X[J]
X2 = (X[J]+X[2]-X[1]):(X[2] - X[1]):(4*X[end])

X = XJ
ho = wave_effs[1].basis_order
S = 1;
species = [specie]
J = length(X) - 1
S_mat = OffsetArray{Complex{T}}(undef, 0:J, -2ho:2ho);
for j = 0:J, m = -2ho:2ho
    S_mat[j,m] = integrate_S(m, X[j+1]; θin = θin)
end

Z = OffsetArray{Complex{T}}(undef,-ho:ho);
for m = 0:ho
    Z[m] = Zn(ω,specie,medium,m)
    Z[-m] = Z[m]
end

Es = [
    [
        sum(
            species[s].num_density * Z[n] * im^T(n+1) * S_mat[J-l,n-m] *
            exp(im*X[J+1]*wave_effs[p].k_eff*cos(wave_effs[p].θ_eff) - im*n*wave_effs[p].θ_eff) *
             wave_effs[p].amplitudes[n+ho+1] / (wave_effs[p].k_eff*cos(wave_effs[p].θ_eff) + cos(θin))
        for n = -ho:ho, s = 1:S)
    for l = 0:J, p in eachindex(wave_effs)]
for m = -ho:ho]
E_mat = vcat(Es...);

σ =  trapezoidal_scheme(X2) # integration scheme: trapezoidal

scalar = 2*specie.num_density/cos(θin)
EE = scalar*[
    sum(
        (-1)^n*im^T(-m)*exp(im*(m-n)*θin)*exp(-im*n*wave_effs[p].θ_eff)*Z[n]*
        wave_effs[p].amplitudes[n+ho+1,1] * exp(-im*xl*cos(θin)) *
        sum(σ.*exp.(im.*X2.*wave_effs[p].k_eff*cos(wave_effs[p].θ_eff)/k .+ im.*X2.*cos(θin)))
    for n = -ho:ho, p in eachindex(wave_effs))
for xl in XJ, m = -ho:ho]
EE2 = sum(E_mat,2)[:]
@test maximum(abs.(EE[:]./EE2 .-1)) < 1.04e-2
