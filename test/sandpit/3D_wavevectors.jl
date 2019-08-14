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
# R_kps = wavenumbers(ω, medium, species; dim = 3, hankel_order = hankel_order,
#     num_wavenumbers=27, tol = 1e-8, mesh_points = 10, mesh_size = 2.0)
R_kps = Complex{Float64}[2.30521+1.85168im, -1.46227+2.03392im, -1.45281+2.03493im, 10.9431+5.60625im, -10.0589+5.61005im, -10.0562+5.61218im, 19.1007+6.88548im, -18.2139+6.88838im, 27.0945+7.7147im, 35.0209+8.32989im, 35.0222+8.33306im, 42.9157+8.82717im, 50.7889+9.23896im, 58.6476+9.58905im, 58.649+9.59209im, 66.5001+9.90123im, 74.3447+10.1762im, 82.1845+10.4237im, 90.0205+10.6489im]

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
