using EffectiveWaves
using Statistics

include("optimise_wavenumber.jl")

## Choose objective function
function f_slow(sps, medium,ωs,kTs)
  mean(ωs./real(kTs)) # average sound speed over whole frequency range
end
function f_constant_speed(sps, medium,ωs,kTs)
  abs(mean(ωs./real(kTs)) - 1000.0) # difference between average sound speed and 1000.0 sound speed
end
f_obj = f_constant_speed # choose which function to minimise

medium = Medium(ρ=1.0,c=1400.0) # choose the background medium
ωs = 0.01:0.01:1.0 # angular frequency range

## How many species to consider
num_species=2

## Choose which fields to optimise
opt_fields = [(:r,(0.001,2.0)),(:volfrac,(0.,0.12))] # this example will vary the radius, :r, in the interval (0.001,2.0), and vary the volume fraction, :volfrac, in the interval (0.,0.12).

# and which fields should take fixed values (should be same lenght as num_species)
fix_fields = [(:ρ,[0.1,100.]),(:c,[2000.0+0.0im,100.0+0.0im])] # this example will use density,:ρ , 0.1 for the first species and 100. for the second species. The analogous applies to the sound speed, :c.

# find the optimal values for the opt_fields.
species = optimal_species(f_obj, medium, ωs;
              opt_fields = opt_fields,
              fix_fields = fix_fields,
              num_species=num_species, MaxTime=100., method = :xnes)
# To choose other methods see variable: opt_methods. To auto select a method, just omit the method option, or equivalently use method = :choose_method.

kTs_arr =[
  wavenumber_low_volfrac(ωs, medium, sps)
for sps in [[species[1]], [species[2]], species] ]

speed_arr = [ ωs./real(kTs) for kTs in kTs_arr]
atten_arr = imag(kTs_arr)

using Plots
unicodeplots()

labs = ["void" "stone" "mix" "medium"];

p1 = plot(ωs./real(medium.c),speed_arr, xlabel="k", ylabel="wave speed", labels=labs,);
p2 = plot(ωs./real(medium.c),atten_arr, xlabel="k", ylabel="attenuation", labels=labs,);
plot(p1,p2)
