using EffectiveWaves
using Plots
pyplot()

medium = Medium(ρ=1.0, c=1.0)
species = [Specie(ρ=0.,r=0.6, c=0.2, volfrac=0.3)]

ωs = [1.0]

ws = wavenumbers(ωs[1], medium, species; dim = 3, num_wavenumbers=10)


wavesystem1 = wavematrix3DPlane(ωs[1], medium, species; hankel_order=4, θ_inc=2.0)

det(wavesystem1(1.0+0.2im))

wavesystem2 = wavematrix3DPlane(ωs[1], medium, species; hankel_order=4, θ_inc=0.0, φ_inc=0.4)
det(wavesystem2(1.0+0.2im))

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
