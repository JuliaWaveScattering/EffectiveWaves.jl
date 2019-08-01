using EffectiveWaves
using LinearAlgebra
using Plots
pyplot()

medium = Medium(ρ=1.0, c=1.0)
species = [Specie(ρ=0.,r=0.6, c=0.2, volfrac=0.3)]

ωs = [1.0]

kp = 1.0+0.4im
kps = collect(0.2:0.2:1.0) .+ collect(0.3:0.2:1.1) .* im
hankel_order = 2

wavesystemPlane = wavematrix3DPlane(ωs[1], medium, species; hankel_order=hankel_order)
detsplanes = det.(wavesystemPlane.(kps))
detP(k) = det(wavesystemPlane(k))

# dets1 = det.(wavesystem1.(kps))
# wavesystem2 = wavematrix3DPlane(ωs[1], medium, species; hankel_order=hankel_order, θ_inc=2.0, φ_inc=0.4)
# det(wavesystem2(kp))

wavesystem3D = wavematrix3D(ωs[1], medium, species; tol = 1e-4, hankel_order=hankel_order);
detR(k) = det(wavesystem3D(k))
dets3D = det.(wavesystem3D.(kps))

ω = 1.
k = ω/medium.c

k0 = real(k)
x = k0 .* LinRange(-5.001,5,150)
y = k0 .* LinRange(0.001,5.0,150)

X = repeat(x',length(y),1)
Y = repeat(y,1,length(x))

ZPs = map( (x,y) -> abs(detP(x+y*im)),X,Y)
ZRs = map( (x,y) -> abs(detR(x+y*im)),X,Y)

h1 = heatmap(x,y,ZPs, xlab = "Re k*", ylab = "Im k*", title="det(I + G) - moderate ω and dirichlet", clims = (0.,0.3))
h2 = heatmap(x,y,ZRs, xlab = "Re k*", ylab = "Im k*", title="det(I + G) - moderate ω and dirichlet", clims = (0.,0.3))



abs.(detsplanes)
abs.(dets3D)

kps = wavenumbers(ωs[1], medium, species; dim = :plane, num_wavenumbers=5)
det.(wavesystemPlane.(kps))

kps3D = wavenumbers(ωs[1], medium, species; dim = 3, num_wavenumbers=5)
det.(wavesystem3D.(kps))
det.(wavesystem3D.(kps3D))
abs.(det.(wavesystemPlane.(kps3D)))
abs.(det.(wavesystemPlane.(kps)))

eff_medium = effective_medium(medium, species)
k_eff_lows = ωs./eff_medium.c

k_eff_φs = wavenumber_low_volfrac(ωs, medium, species)
# num_wavenumbers =1 almost always finds the wavenubmer with the smallest attenuation

k_effs_arr = [
    wavenumbers(ω, medium, species; tol=tol, num_wavenumbers=1)
for ω in ωs]

inds = [argmin(abs.(k_effs_arr[i] .- k_eff_φs[i])) for i in eachindex(ωs)]
k_effs2 = [k_effs_arr[i][inds[i]] for i in eachindex(inds)]


" Returns all the transmitted effective wavenumbers"
wavenumbers(ω::T, medium::Medium{T}, specie::Specie{T}; kws...) where T<:Number = wavenumbers(ω, medium, [specie]; kws...)
