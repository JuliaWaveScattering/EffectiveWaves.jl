using EffectiveWaves
using LinearAlgebra
using Plots
pyplot()

medium = Medium(ρ=1.0, c=1.0)
species = [Specie(ρ=0.1,r=0.4, c=0.2, volfrac=0.3)]

ωs = [1.0]
ω = ωs[1]
k = ω/medium.c
hankel_order = ho = 2


kps = wavenumbers(ω, medium, species; symmetry = :plane, dim = 3,
    hankel_order = hankel_order, num_wavenumbers=17,
    tol = 1e-8,
    mesh_points = 14, mesh_size = 2.0)

# Doesnt really work
R_kps = wavenumbers(ω, medium, species; dim = 3, hankel_order = hankel_order,
    num_wavenumbers=27, tol = 1e-8, mesh_points = 10, mesh_size = 2.0)


wavesystem2D = wavematrix2D(ωs[1], medium, species; hankel_order=hankel_order);
det2D(kp) = det(wavesystem2D(kp))
kps2 = wavenumbers(ω, medium, species; dim = 2, hankel_order = hankel_order,
    num_wavenumbers=27, tol = 1e-8, mesh_points = 10, mesh_size = 2.0);


wavesystemPlane = wavematrix3DPlane(ωs[1], medium, species; hankel_order=hankel_order)
detP(kp) = det(wavesystemPlane(kp))
detP.(kps)

inds = findall( abs.(detP.(kps)) .< 1e-5);
kps = kps[inds]

wavesystem3D = wavematrix3D(ωs[1], medium, species; hankel_order=hankel_order);
detR(kp) = (det(wavesystem3D(kp)))^(1/2)
detR.(kps)


MRsvd = svd(wavesystem3D(kps[2]))
MRsvd = svd(wavesystem3D(kps[1]))
MRsvd.S
wave_vec =  MRsvd.V[:,end]
wave_mat = reshape(wave_vec,((ho+1)^2,(ho+1)^2))
