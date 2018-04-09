include("low_volfrac.jl")


#  Specie(ρ,r; volfrac=0.1), where ρ = density and r = radius.
#  Specie(ρ,r,c,num_den)
# sps = [ Specie(ρ=2200.,r=30.e-6,c=8433.,volfrac=0.1), Specie(ρ=Inf,r=100.0e-6,c=2.0,volfrac=0.2)]
sp = Specie(ρ=2200.,r=30.e-6,c=8433.,volfrac = 0.1)
sp = Specie(ρ=2200,r=30e-6,c=8433.,volfrac = 0.1)

ω = 30.0e6
# 2pi*2250.0/(ω*30e-6)

medium = Medium(ρ=2329.0,c = 2250.0)
kT = sqrt(wavenumber_low_volfrac(ω, medium, sp; verbose=true, rel_tol=0.00002))

ωs = linspace(0.01,1.0,60)*30.0e6
kTs = sqrt(wavenumber_low_volfrac(ωs, medium, sp))

speed_arr = ωs./real(kTs)
atten_arr = imag(kTs)

using Plots
pyplot()

p1 = plot(ωs./real(medium.c),speed_arr,  ylabel="wave speed", xlabel="k");
p2 = plot(ωs./real(medium.c), atten_arr, ylabel="attenuation", xlabel="k");
plot(p1,p2)

## An example where we vary the species

height=450
pyplot(linewidth=3, size=(2*height,height))
# pyplot(linewidth=3)
Plots.scalefontsizes(1.5)

# for fixed total volfrac fraction
medium = Medium(ρ=1.0,c = 1.0)
ωs = linspace(0.01,1.0,60)
volfrac = 0.25
r1 = 0.5
r2 = 1.5
N=5
vols = linspace(0.0,1.0,N)*volfrac

kTs_arr = [
  begin
    sp1 = Specie(0.0, r1; volfrac=vols[i])
    sp2 = Specie(Inf, r1; volfrac=volfrac-vols[i])
    [ sqrt(wavenumber_low_volfrac(ω, medium, [sp1,sp2])) for ω in ωs]
  end
for i = 1:N];

speed_arr = [ ωs./real(kTs) for kTs in kTs_arr]
atten_arr = imag(kTs_arr)

labs = reshape( map(v -> "void vol = $(Int(round(100*v)))%",vols),1, length(vols));
p1 = plot(ωs./real(medium.c), speed_arr, labels=labs, ylabel="wave speed" ,xlabel="k");
p2 = plot(ωs./real(medium.c), atten_arr, labels=labs, xlabel="k", ylabel="attenuation");

plot(p1,p2)
savefig("varying_voids.png")

Plots.scalefontsizes(1./1.5)
