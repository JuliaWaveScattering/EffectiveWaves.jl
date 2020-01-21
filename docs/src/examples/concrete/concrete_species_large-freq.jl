# This example is used to generate plots in the paper, "Reflection from a multi-species material and its transmitted effective wavenumber." Proc. R. Soc. (2018): 20170864.

# Everything related to ploting has been commented so that this example can run with requiring the Plots package. That way, this example can be checked everytime the package is modified.

using EffectiveWaves

# using LaTeXStrings
# using Plots
# height=500
 # pyplot(linewidth=3, size=(2.6*height,height), border=false)

 # Plots.scalefontsizes(1.7)

## choose material

  # concrete
  medium = LimeStone
  inclusion1 = AirDry
  inclusion2 = Brick

  ωfactor = 50.0;
  ωs = 20*ωfactor.*LinRange(real(medium.c/10000),real(medium.c),400) # k from 0 to 1

  volfrac = 0.16
  r1 = 0.1/ωfactor; vol1 = 0.06
  r2 = 1.0/ωfactor; vol2 = volfrac - vol1

  sp1 = Specie(ρ=inclusion1.ρ ,r=r1, c=inclusion1.c, volfrac = vol1)
  sp2 = Specie(ρ=inclusion2.ρ ,r=r2, c=inclusion2.c, volfrac = volfrac-vol1)

  # True wavenumber
  kTs  = wavenumber_low_volumefraction(ωs, medium, [sp1,sp2]; tol=0.5e-5);
  # Approximate challis
  kTCs = wavenumber_challis(ωs, medium, [sp1,sp2]; basis_order=10);

  kTs_arr = [kTs,kTCs];
  speed_arr = [ ωs./real(ks) for ks in kTs_arr];
  push!(speed_arr, 0 .*ωs .+ real(medium.c))
  atten_arr = imag.(kTs_arr)

  styles = [:solid :dashdot :dot]
  # labs = [L"Effective $k_{*}$" L"Approximate $k_{*C}$" "Lime Stone"]
  ys_arr = speed_arr;
  xs = r1.*(ωs./real(medium.c));
  m =5;
  y1 = min(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  y2 = max(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  # p1 = plot(xs, ys_arr, xlabel=L"a_S k", ylabel="sound speed (m/s)", labels=labs, line = styles)

  ys_arr = atten_arr;
  # labs = [L"Effective $k_{*}$" L"Approximate $k_{*C}$"]
  # p2 = plot(xs, ys_arr, labels=labs, xlabel=L"a_S k", ylabel="attenuation (1/m)", line=styles)

  # plot(p1,p2)
  # try mkdir("media") end
  # savefig("media/compare_concrete_large-w.png")
  # savefig("media/compare_concrete_large-w.pdf")
  # gui()

# Plots.scalefontsizes(1/1.7)
