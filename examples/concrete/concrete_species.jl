using EffectiveWaves
using LaTeXStrings

## choose material

  # concrete
  medium = LimeStone
  inclusion1 = AirDry
  inclusion2 = Brick

  ωfactor = 50.0;
  ωs = ωfactor.*linspace(real(medium.c/10000),real(medium.c),100) # k from 0 to 1

  volfrac = 0.16
  r1 = 0.1/ωfactor; vol1 = 0.06
  r2 = 1.0/ωfactor; vol2 = volfrac - vol1

  sp1 = Specie(ρ=inclusion1.ρ ,r=r1, c=inclusion1.c, volfrac = vol1)
  sp2 = Specie(ρ=inclusion2.ρ ,r=r2, c=inclusion2.c, volfrac = volfrac-vol1)

  # True wavenumber
  kTs  = wavenumber_low_volfrac(ωs, medium, [sp1,sp2]);
  # Approximate wavenumber
  kTLSs = two_species_approx_wavenumber(ωs, medium, [sp1,sp2]);
  # Approximate challis
  kTCs = wavenumber_challis(ωs, medium, [sp1,sp2]);

  kTs_arr = [kTs,kTLSs,kTCs];
  speed_arr = [ ωs./real(ks) for ks in kTs_arr];
  push!(speed_arr, 0.*ωs + real(medium.c))
  atten_arr = imag.(kTs_arr)

  styles = [:solid :dashdot :dashdot :dot]
  labs = [L"Effective $k_{*}$" L"Approximate $k_{*LS}$" L"Approximate $k_{*C}$" "Lime Stone"]
  ys_arr = speed_arr;
  xs = r1.*(ωs./real(medium.c));
  m =5;
  y1 = min(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  y2 = max(ys_arr[1][1:m]..., ys_arr[2][1:m]...);

  using Plots
  height=500
  # unicodeplots()
  pyplot(linewidth=3, size=(2.6*height,height), border=false)

  Plots.scalefontsizes(1.7)
  plot(xs, ys_arr, xlabel=L"a_S k", ylabel="sound speed (m/s)", labels=labs, line = styles, xlims = (-0.002,maximum(xs))
        , ylims = ( min(ys_arr[1]...,ys_arr[3]...,ys_arr[4]...)*0.995,  max(ys_arr[1]...,ys_arr[3]...,ys_arr[4]...)*1.005));
  p1 = gray_square([xs[1],xs[m]],[y1,y2],plot!);

  ys_arr = atten_arr;
  labs = [L"Effective $k_{*}$" L"Approximate $k_{*LS}$" L"Approximate $k_{*C}$"]
  plot(xs, ys_arr, labels=labs, xlabel=L"a_S k", ylabel="attenuation (1/m)", line=styles, xlims = (-0.002,maximum(xs))
        , ylims = ( min(ys_arr[1]...,ys_arr[3]...)*0.995,  max(ys_arr[1]...,ys_arr[3]...)*1.005));
  y1 = min(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  y2 = max(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  p2 = gray_square([xs[1],xs[m]],[y1,y2],plot!);

  plot(p1,p2)
  try mkdir("media") end
  savefig("media/compare_concrete.png")
  savefig("media/compare_concrete.pdf")
  gui()

## Zoomed in version
  ωs = linspace(ωs[1],ωs[m],250) # k from 0 to 1
  m = length(ωs);
  kTs  = wavenumber_low_volfrac(ωs, medium, [sp1,sp2]);
  kTLSs = two_species_approx_wavenumber(ωs, medium, [sp1,sp2]);
  kTCs  = wavenumber_challis(ωs, medium, [sp1,sp2]);

  kTs_arr = [kTs,kTLSs,kTCs];
  speed_arr = [ ωs./real(ks) for ks in kTs_arr];
  atten_arr = imag.(kTs_arr)

  labs = [L"Effective $k_{*}$" L"Approximate $k_{*LS}$" L"Approximate $k_{*C}$"]

  ys_arr = speed_arr;
  xs = r1.*(ωs./real(medium.c));
  plot(xs, ys_arr, xlabel=L"a_S k", ylabel="sound speed (m/s)", labels=labs
                 , border = false, line = styles, xlims = (0,maximum(xs)));
  y1 = min(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  y2 = max(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  p1 = gray_square([xs[1],xs[m]],[y1,y2],plot!);

  ys_arr = atten_arr;
  labs = [L"Effective $k_{*}$" L"Approximate $k_{*LS}$" L"Approximate $k_{*C}$"]
  plot(xs, ys_arr, xlabel=L"a_S k", ylabel="attenuation (1/m)", labels=labs
                 , border=false, line = styles, xlims = (0,maximum(xs)));
  y1 = min(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  y2 = max(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  dy = abs(y2- y1)/120.0;
  p2 = gray_square([xs[1],xs[m]],[y1,y2],plot!);

  plot(p1,p2)
  savefig("media/compare_concrete_zoom.png")
  savefig("media/compare_concrete_zoom.pdf")

Plots.scalefontsizes(1/1.7)
