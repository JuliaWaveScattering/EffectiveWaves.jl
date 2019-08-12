using EffectiveWaves
using LinearAlgebra
using Plots
pyplot()

medium = Medium(ρ=1.0, c=1.0)
species = [Specie(ρ=0.1,r=0.01, c=0.2, volfrac=0.3)]

ωs = [0.001]
ω = ωs[1]
k = ω/medium.c


kp = 1.0+0.4im
kps = collect(0.2:0.2:1.0) .+ collect(0.3:0.2:1.1) .* im
hankel_order = ho = 2

wavesystemPlane = wavematrix3DPlane(ωs[1], medium, species; hankel_order=hankel_order)
detP(kp) = det(wavesystemPlane(kp))

# detsplanes = det.(wavesystemPlane.(kps))
# MPs = wavesystemPlane.(kps)

wavesystem3D = wavematrix3D(ωs[1], medium, species; tol = 1e-4, hankel_order=hankel_order);
detR(kp) = det(wavesystem3D(kp))

# wavesystem3D_2 = wavematrix3D_allocate(ωs[1], medium, species; tol = 1e-4, hankel_order=hankel_order);
# detR2(kp) = det(wavesystem3D_2(kp))

# @time detR.(kps);
#   1.996288 seconds (469.46 k allocations: 52.257 MiB, 1.06% gc time)
#
# @time detR.(kps);
#   2.103966 seconds (327.72 k allocations: 45.012 MiB, 0.39% gc time)
#
# @time detR.(kps);
# 0.084732 seconds (32.86 k allocations: 3.510 MiB)


# Effective medium
    dim = 3; T= Float64
    β = medium.ρ*medium.c^2 # medium bulk modulus
    φ = sum(volume_fraction.(species; dim = dim))
    β_eff = one(T)/((1-φ)/β + sum(volume_fraction(s; dim = dim)/(s.ρ*s.c^2) for s in species))
    ρ_frac = sum(volume_fraction(s; dim = dim)*(medium.ρ - s.ρ)/(medium.ρ + 2*s.ρ) for s in species)

    ρ_eff = medium.ρ*(1 - ρ_frac)/(1 + T(dim - 1) * ρ_frac)
    eff_medium = Medium(ρ=ρ_eff, c=sqrt(β_eff/ρ_eff))
    k_low1 = ω/eff_medium.c

    ρ_eff = medium.ρ*(1 - ρ_frac)/(1 + (ho/(2ho+3))*ρ_frac)
    eff_medium = Medium(ρ=ρ_eff, c=sqrt(β/ρ_eff))
    k_low2 = ω/eff_medium.c

    ρ_eff = medium.ρ*(1 - ρ_frac)/(1 + ((ho-1)/(2ho+1))*ρ_frac)
    eff_medium = Medium(ρ=ρ_eff, c=sqrt(β_eff/ρ_eff))
    k_low3 = ω/eff_medium.c

kps = wavenumbers(ω, medium, species; dim = 3,
    num_wavenumbers=7, tol = 1e-8, mesh_points = 10, mesh_size = 2.0)
detP.(kps)
detR.(kps)

k_lows = [k_low1,k_low2,k_low3]
k1 = minimum(real.(k_lows[[1,3]]))
k2 = maximum(real.(k_lows))
x = LinRange(0.9*k1,1.1*k2,110)
y = k2 .* LinRange(-0.01,0.2,90)
# y = LinRange(-0.0015,0.0015,90)

X = repeat(x',length(y),1)
Y = repeat(y,1,length(x))

ZPs = map( (x,y) -> abs(detP(x+y*im)),X,Y)
ZRs = map( (x,y) -> abs(detR(x+y*im)),X,Y)

minR = 1e-18
minP = 0.2

ZRs2 = [ (abs(ZRs[ind]) < minR ? ZRs[ind] : NaN ) for ind in LinearIndices(ZRs)];
ZPs2 = [ (abs(ZPs[ind]) < minP ? ZPs[ind] : NaN ) for ind in LinearIndices(ZPs)];

# h1 = heatmap(x,y,ZPs, xlab = "Re k*", ylab = "Im k*", title="Planewave det(I + G), ka~1", clims = (0.,1.0))

h2 = heatmap(x,y,ZRs2, xlab = "Re k*", ylab = "Im k*", title="General det(I + G), ka~$(ω*species[1].r)", clims = (0.,minR))
h1 = heatmap(x,y,ZPs2, xlab = "Re k*", ylab = "Im k*", title="Planewave det(I + G), ka~$(ω*species[1].r)", clims = (0.,minP))

# Z = map( (x,y) -> abs(((x^2-1)^2*(y*im - 1)*(y*im - 2))),X,Y)
# Z = [ (abs(Z[ind]) < 1.0 ? Z[ind] : NaN ) for ind in LinearIndices(Z)]
# heatmap(x,y,Z, xlab = "Re k*", ylab = "Im k*", clims = (0.,1.0))

abs.(detsplanes)
abs.(dets3D)

kps = wavenumbers(ωs[1], medium, species; symmetry = :plane, dim = 3, num_wavenumbers=5)
detP.(kps)

MPsvd = svd(wavesystemPlane(kps[1]))
MPsvd.S
MPsvd.V[:,end]

MRsvd = svd(wavesystem3D(kps[1]))
MRsvd = svd(wavesystem3D(kps[2]))
MRsvd.S
MRsvd.V[:,end]

v = rand(3)

M = [v v v]
Msvd = svd(M)
Msvd.S
Msvd.V[:,3]

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
