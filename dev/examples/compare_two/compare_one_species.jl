using EffectiveWaves

using Plots
# unicodeplots() # alternative plot backend
height=500
pyplot(linewidth=3, size=(2.0*height,height))

# for fixed total volfrac fraction
medium = Medium(ρ=3200.0,c = 2595.0)
ωs = 2.0*pi*LinRange(1.0e1,1.0e7,100)

volfrac = 0.12
r1 = 30.0e-8
r2 = 2*30.0e-6

sp1 = Specie(ρ=2329.,r=r1,c=8433.,volfrac = 0.06)
sp2 = Specie(ρ=2329.,r=r2,c=8433.,volfrac = volfrac-0.06)

vol = pi*sp1.num_density*sp1.r^2

rhoeff =  medium.ρ*(medium.ρ + sp1.ρ - vol*(medium.ρ - sp1.ρ))/
(medium.ρ + sp1.ρ + vol*(medium.ρ - sp1.ρ))

# rhoeff =  medium.ρ*(1.0-vol) + vol*sp1.ρ
# rhoeff =  medium.ρ*(1 + volfrac)/(1 - volfrac)

kTs = [ wavenumber_low_volumefraction(ω, medium, [sp1,sp2]) for ω in ωs]
kTaprx =[ two_species_approx_wavenumber(ω, medium, [sp1,sp2]) for ω in ωs]

kTSs =[ wavenumber_low_volumefraction(ω, medium, [sp1]) for ω in ωs]
kTLs = [
  begin
    mS = Medium(ρ=rhoeff, c= ωs[i]/kTSs[i])
    wavenumber_low_volumefraction(ωs[i], mS, [sp2])
  end
for i in eachindex(ωs)];

# speed_arr = [ [ωs./real(kTs)], [0.029 .+ ωs./real(kTLs)]]
speed_arr = [ ωs./real(kTs), ωs./real(kTLs), ωs./real(kTaprx)]
atten_arr = imag([kTs,kTLs, kTaprx])

labs = reshape(["kT","kTL", "kTaprx"] ,1, 3);
p1 = plot(r1.*(ωs./real(medium.c))/(2pi), speed_arr, labels=labs, xlabel="a/λ", ylabel="wave speed");
p2 = plot(r1.*(ωs./real(medium.c))/(2pi), atten_arr, labels=labs, xlabel="a/λ", ylabel="attenuation");
plot(p1,p2)
gui()

speed_arr = [abs.(1 .- real(kTLs)./real(kTs))]
p1 = plot(r1.*(ωs./real(medium.c))/(2pi), speed_arr, labels=labs, xlabel="a/λ", ylabel="relative wave speed");
