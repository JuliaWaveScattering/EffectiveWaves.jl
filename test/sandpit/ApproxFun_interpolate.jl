using ApproxFun
using EffectiveWaves

T = Float64
θin = 0.13
k=1.; ho = 2
medium = Medium(1.0,1.0+0.0im)
ω = real(k*medium.c)
specie = Specie(ρ=0.1,r=0.1, c=0.5, volfrac=0.1)

# Guess a course mesh
k_eff0 = wavenumber_low_volfrac(ω, medium, [specie]; tol = 1e-12)
max_x = 10 .* k/imag(k_eff0)
x = 0.0:0.1:max_x

wave0 = EffectiveWave(ω, k_eff0, medium, [specie]; θin = θin, hankel_order = ho, tol=1e-8)
wave_avg0 = AverageWave(x, wave0)

R_eff = reflection_coefficient(ω, wave0, medium, [specie]; θin = θin, hankel_order = ho)
R = reflection_coefficient(ω, wave_avg0, medium, specie; θin = θin)

amps = wave_avg0;
# k = ω/medium.c
M = amps.hankel_order
σ = trap_scheme(amps.x);

R_ms = T(2)*specie.num_density/(cos(θin)*k^2)*
[
    sum(
        im^T(m)*exp(-im*θin*m)*amps.amplitudes[j,m+M+1,1]*exp(im*amps.x[j]*cos(θin))*σ[j]
    for j in eachindex(amps.x))
for m = -M:M];
sum(R_ms) - R

x = 0.:0.02:10

# function space
S = Chebyshev(minimum(x)..maximum(x));
n1 = length(x); m1 = Int(round(n*0.4));
V = Array{Float64}(n1,m1); # Vandermonde matrix
for k = 1:m1
    V[:,k] = Fun(S,[zeros(k-1);1]).(x)
end

v = exp(x.*(-1.0 + 1.0im))
σ = trap_scheme(x)
sum(σ[i]*v[i] for i in eachindex(x))

f = Fun(S,V\v)
sum(Fun(S,V\v))
1/(1.0+1.0im)

R_funs = map(-M:M) do m
    v = [
        im^T(m)*exp(-im*θin*m)*amps.amplitudes[j,m+M+1,1]*exp(im*amps.x[j]*cos(θin))
    for j in eachindex(amps.x)]
    sum(Fun(S,V\v))
end

R_fun = sum(R_funs)*T(2)*specie.num_density/(cos(θin)*k^2)


    R_m = sum(Fun(S,V\v));
sum(R_m)
