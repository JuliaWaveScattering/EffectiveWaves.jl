include("../../src/EffectiveWaves.jl")

using LaTeXStrings
using Plots
height=500
# unicodeplots()
 pyplot(linewidth=3, size=(2.6*height,height), border=false)

 Plots.scalefontsizes(1.8)

filename="fluid"
mediumname = "water"
## choose material

  # fluid
  medium = WaterDistilled
  # inclusion1 = SodiumNitrate # surfactant ?
  inclusion1 = Glycerol # surfactant ?
  inclusion2 = Hexadecane

  # medium = Hexadecane
  # # inclusion1 = SodiumNitrate # surfactant ?
  # inclusion1 = SodiumNitrate # surfactant ?
  # inclusion2 = Glycerol

  ωfactor = 4000;
  ωs = 20*ωfactor.*linspace(real(medium.c/10000),real(medium.c),200) # k from 0 to 1

  volfrac = 0.22
  r1 = 0.1/ωfactor; vol1 = 0.11;
  r2 = 1.0/ωfactor; vol2 = volfrac - vol1

  sp1 = Specie(ρ=inclusion1.ρ, r=r1, c=inclusion1.c, volfrac = vol1)
  sp2 = Specie(ρ=inclusion2.ρ, r=r2, c=inclusion2.c, volfrac = volfrac-vol1)

  # True wavenumber
  kTs  = multispecies_wavenumber(ωs, medium, [sp1,sp2]);

  # Approximate Challis wavenumber
  kTCs = multispecies_challis(ωs, medium, [sp1,sp2]; hankel_order=5);

  speed_arr = [ ωs./real(kTs), ωs./real(kTCs), 0.*ωs + real(medium.c)]
  atten_arr = imag([kTs,kTCs])

  styles = [:solid :dashdot :dot]
  labs = [L"Effective $k_{*}$" L"Approximate $k_{*C}$" mediumname];
  ys_arr = speed_arr;
  xs = r1.*ωs./real(medium.c); xlabel = L"k a_S";
  m =5;
  p1= plot(xs, ys_arr, xlabel=xlabel, ylabel="sound speed (m/s)", labels=labs, line = styles
            , ylims = (minimum(minimum.(ys_arr))*0.995, maximum(maximum.(ys_arr))*1.005));

  ys_arr = atten_arr;
  labs = [L"Effective $k_{*}$" L"Approximate $k_{*C}$"];
  p2 = plot(xs, ys_arr, labels=labs, xlabel=xlabel, ylabel="attenuation (1/m)", line = styles
              , ylims = (minimum(minimum.(ys_arr))*0.995, maximum(maximum.(ys_arr))*1.005));
  plot(p1,p2)
  savefig("compare_$(filename)_large-w.png")
  savefig("compare_$(filename)_large-w.pdf")

  Plots.scalefontsizes(1/1.8)
