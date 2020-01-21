using EffectiveWaves

include("../discrete_wave/integral_form.jl")

# physical parameters
medium = Medium(1.0,1.0+0.0im)
specie = Specie(ρ=0.5, r=0.5, c=0.2, volfrac=0.1)

k = 0.8;
θin = 0.0
ω = real(k*medium.c)

ho = maximum_basis_order(ω, medium, [specie]; tol=tol*1e4)

X0 = 6.
X1 = 8.
amps = map(8:1:16) do X_max
    X = 0.0:0.05:X_max
    i_X0 = findmin(abs.(X0 .- X))[2]
    i_X1 = findmin(abs.(X1 .- X))[2]
    avg_wave = DiscretePlaneWaveMode(ω, medium, specie; basis_order=ho, X=X)
    avg_wave.amplitudes[i_X0:i_X1,:,:]
end

errs = [norm(amps[i][:] - amps[end][:])/norm(amps[end][:]) for i in eachindex(amps)[1:end-1]]
sum(log.(errs) .< -2:-1.2:-10.5) == 8

# plot((8:1:16)[1:end-1],errs)
