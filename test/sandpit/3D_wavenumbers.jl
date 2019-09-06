using EffectiveWaves
using LinearAlgebra
using Plots
pyplot()

medium = Medium(ρ=1.0, c=1.0)
species = [Specie(ρ=0.1,r=1.0, c=0.01, volfrac=0.3)]

ωs = [0.5]
ω = ωs[1]
k = ω/medium.c
tol = 1e-8

kp = 1.0+0.4im
kps = collect(0.2:0.2:1.0) .+ collect(0.3:0.2:1.1) .* im
hankel_order = ho = 2

wavesystemPlane = wavematrix3DPlane(ωs[1], medium, species; tol = 1e-5, hankel_order=hankel_order)
detP(kp) = det(wavesystemPlane(kp))

# dispersion = dispersion_function(ω, medium, species; tol = low_tol, dim=dim, symmetry = :plane, hankel_order = ho)

# dim = 3
# k_effs = wavenumbers_path(ω, medium, species;
#     num_wavenumbers=7, dim=dim, symmetry = :plane,
#     hankel_order = ho, tol=tol)

# detsplanes = det.(wavesystemPlane.(kps))
# MPs = wavesystemPlane.(kps)

wavesystem3D = wavematrix3D(ωs[1], medium, species; tol = 1e-5, hankel_order=hankel_order);
detR(kp) = sqrt(det(wavesystem3D(kp)))

wavesystem3D = wavematrix3D(ωs[1], medium, species; tol = 1e-5, hankel_order=1);
detR(kp) = sqrt(det(wavesystem3D(kp)))

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

kp0s = wavenumbers(ω, medium, species; dim = dim, symmetry = :plane, hankel_order = ho,
    num_wavenumbers=20, tol = tol, mesh_points = 14, mesh_size = 3.0)

kp_vecs = reduce_kvecs([[real(kp),imag(kp)] for kp in kp0s], 1e-5)
kps = [kp[1]+im*kp[2] for kp in kp_vecs]
kps = kps[1:30]

# inds = findall(abs.(detP.(kps)) .< 1e-6)
# kps = kp0s[inds]

# For low frequency
    # k_lows = [k_low1,k_low2,k_low3]
    # k1 = minimum(real.(k_lows[[1,3]]))
    # k2 = maximum(real.(k_lows))
    # x = LinRange(0.9*k1,1.1*k2,110)
    # y = k2 .* LinRange(-0.01,0.2,90)

# Given some kps
    # dk1 = real(kps[1]) / 4
    # dk2 = imag(kps[1]) / 4
    k1 = maximum(abs.(real.(kps)))
    k2 = imag(kps[end])
    x = LinRange(-k1,k1,6600)
    y = LinRange(-0.01,k2,160)

X = repeat(x',length(y),1)
Y = repeat(y,1,length(x))

ZPs = map( (x,y) -> abs(detP(x+y*im)),X,Y)
# ZRs = map( (x,y) -> abs(detR(x+y*im)),X,Y)

minR = 1e-18
minP = 0.8

# ZRs2 = [ (abs(ZRs[ind]) < minR ? ZRs[ind] : NaN ) for ind in LinearIndices(ZRs)];
ZPs2 = [ (abs(ZPs[ind]) < minP ? ZPs[ind] : NaN ) for ind in LinearIndices(ZPs)];

# h1 = heatmap(x,y,ZPs, xlab = "Re k*", ylab = "Im k*", title="Planewave det(I + G), ka~1", clims = (0.,1.0))

# h2 = heatmap(x,y,ZRs2, xlab = "Re k*", ylab = "Im k*", title="General det(I + G), ka~$(ω*species[1].r)", clims = (0.,minR))
h1 = heatmap(x,y,ZPs2, xlab = "Re k*", ylab = "Im k*", title="Planewave det(I + G), ka~$(ω*species[1].r)", clims = (0.,minP))
scatter!(kps)


# Z = map( (x,y) -> abs(((x^2-1)^2*(y*im - 1)*(y*im - 2))),X,Y)
# Z = [ (abs(Z[ind]) < 1.0 ? Z[ind] : NaN ) for ind in LinearIndices(Z)]
# heatmap(x,y,Z, xlab = "Re k*", ylab = "Im k*", clims = (0.,1.0))

abs.(detsplanes)
abs.(dets3D)

h1 = heatmap(x,y,ZPs, xlab = "Re k*", ylab = "Im k*", title="Planewave det(I + G), ka~$(ω*species[1].r)", clims = (0.,1.0))
scatter!(kps)
