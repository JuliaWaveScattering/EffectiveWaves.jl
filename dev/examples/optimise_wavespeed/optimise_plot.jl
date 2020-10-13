include("low_volumefraction.jl")

using JLD
using Plots
plotly()

strfs = ["f_slow1","f_slow2","f_fast1","f_fast2"];
f = strfs[4]
species = load("data/$f.jld")["species"]
ωs = load("data/$f.jld")["ωs"]
medium = load("data/$f.jld")["medium"]

kTs_arr = [
    [ sqrt(wavenumber_low_volumefraction(ω,medium, sps)) for ω in ωs]
for sps in [[species[1]],[species[2]],species]];

speed_arr = [ ωs./real(kTs) for kTs in kTs_arr]
atten_arr = imag(kTs_arr)

labs = ["void" "stone" "mix"];
p1 = plot(ωs./real(medium.c),speed_arr, xlabel="k", ylabel="wave speed", labels=labs,);
p2 = plot(ωs./real(medium.c),atten_arr, xlabel="k", ylabel="attenuation", labels=labs,);
plot(p1,p2)
