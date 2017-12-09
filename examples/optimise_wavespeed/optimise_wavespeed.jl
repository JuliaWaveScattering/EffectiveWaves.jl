include("optimse_two.jl")

using JLD

medium = Medium{Float64}(1.0,1.)
ωs = 0.01:0.01:1.0
MaxTime = 10000.0

function f_slow1(sps::Array{Specie{Float64}}, medium::Medium{Float64},ωs,kTs)
  N = Int(round(length(ωs)/2))
  mean(imag(kTs[1:N])) + mean(ωs[1:N]./real(kTs[1:N]))
end
function f_slow2(sps::Array{Specie{Float64}}, medium::Medium{Float64},ωs,kTs)
  N = Int(round(length(ωs)/2))
  mean(imag(kTs[N:end])) + mean(ωs[N:end]./real(kTs[N:end]))
end
function f_fast1(sps::Array{Specie{Float64}}, medium::Medium{Float64},ωs,kTs)
  N = Int(round(length(ωs)/2))
  mean(imag(kTs[1:N])) - mean(ωs[1:N]./real(kTs[1:N]))
end
function f_fast2(sps::Array{Specie{Float64}}, medium::Medium{Float64},ωs,kTs)
  N = Int(round(length(ωs)/2))
  mean(imag(kTs[N:end])) - mean(ωs[N:end]./real(kTs[N:end]))
end

fs = [f_slow1,f_slow2,f_fast1,f_fast2]

for f in fs
  species = optimal_species(f,medium,ωs; fix_fields = [(:c,[1.0+0.0im,1.0+0.0im])],
              opt_fields = [(:r,(0.01,2.0)),(:volfrac,(0.,0.14)),(:ρ,(0.0,100.))],
              num_species=2, MaxTime=MaxTime, method=:xnes)
  save("data/$f.jld","objective", string(f), "species", species, "ωs", ωs, "medium", medium)
end

# using Plots
# plotly()
# strfs = ["f_slow1","f_slow2","f_fast1","f_fast2"];
# species = load("data/$f.jld")["species"]
#
# kTs_arr = [
#     [ sqrt(multispecies_wavenumber(ω,medium, sps)) for ω in ωs]
# for sps in [[species[1]],[species[2]],species]];
#
# speed_arr = [ ωs./real(kTs) for kTs in kTs_arr]
# atten_arr = imag(kTs_arr)
#
# labs = ["void" "stone" "mix"];
# p1 = plot(ωs./medium.c,speed_arr, xlabel="k", ylabel="wave speed", labels=labs,);
# p2 = plot(ωs./medium.c,atten_arr, xlabel="k", ylabel="attenuation", labels=labs,);
# plot(p1,p2)
