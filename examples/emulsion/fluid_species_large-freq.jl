# This example is used to generate plots in the paper, "Reflection from a multi-species material and its transmitted effective wavenumber." Proc. R. Soc. (2018): 20170864.

# Everything related to ploting has been commented so that this example can run with requiring the Plots package. That way, this example can be checked everytime the package is modified.
using EffectiveWaves

# using LaTeXStrings
# using Plots
# height=500
 # pyplot(linewidth=3, size=(2.6*height,height), border=false)

 # Plots.scalefontsizes(1.8)

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
  kTs  = wavenumber_low_volfrac(ωs, medium, [sp1,sp2]);

  # Approximate Challis wavenumber
  kTCs = wavenumber_challis(ωs, medium, [sp1,sp2]; hankel_order=5);

  speed_arr = [ ωs./real(kTs), ωs./real(kTCs), 0.*ωs + real(medium.c)]
  atten_arr = imag([kTs,kTCs])

  styles = [:solid :dashdot :dot]
  # labs = [L"Effective $k_{*}$" L"Approximate $k_{*C}$" mediumname];
  # xlabel = L"k a_S";
  ys_arr = speed_arr;
  xs = r1.*ωs./real(medium.c);
  # m =5;
  # p1= plot(xs, ys_arr, xlabel=xlabel, ylabel="sound speed (m/s)", labels=labs, line = styles
  #           , ylims = (minimum(minimum.(ys_arr))*0.995, maximum(maximum.(ys_arr))*1.005));

  ys_arr = atten_arr;
  # labs = [L"Effective $k_{*}$" L"Approximate $k_{*C}$"];
  # p2 = plot(xs, ys_arr, labels=labs, xlabel=xlabel, ylabel="attenuation (1/m)", line = styles
  #             , ylims = (minimum(minimum.(ys_arr))*0.995, maximum(maximum.(ys_arr))*1.005));
  # plot(p1,p2)
  # try mkdir("media") end
  # savefig("media/compare_$(filename)_large-w.png")
  # savefig("media/compare_$(filename)_large-w.pdf")
  #
  # Plots.scalefontsizes(1/1.8)
