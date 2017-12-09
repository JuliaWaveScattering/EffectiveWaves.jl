include("multi-species.jl")
include("multi-species_challis.jl")
include("two_species_approximate.jl")
include("materials.jl")
include("graphics.jl")

using LaTeXStrings
using Plots
height=500
# unicodeplots()
 pyplot(linewidth=3, size=(2.6*height,height), border=false)

 Plots.scalefontsizes(1.7)

## choose material

  # concrete
  medium = LimeStone
  inclusion1 = AirDry
  inclusion2 = Brick

  ωfactor = 50.0;
  ωs = ωfactor.*linspace(real(medium.c/10000),real(medium.c),100) # k from 0 to 1
  # ωs = linspace(real(medium.c/10000),real(medium.c/1000),100) # k from 0 to 1
  # ωs = 2.0*pi*linspace(1.0e1,1.0e7,100)

  volfrac = 0.16
  r1 = 0.1/ωfactor; vol1 = 0.06
  r2 = 1.0/ωfactor; vol2 = volfrac - vol1

  sp1 = Specie(ρ=inclusion1.ρ ,r=r1, c=inclusion1.c, volfrac = vol1)
  sp2 = Specie(ρ=inclusion2.ρ ,r=r2, c=inclusion2.c, volfrac = volfrac-vol1)

  # True wavenumber
  kTs  = multispecies_wavenumber(ωs, medium, [sp1,sp2]);
  # Approximate wavenumber
  kTLSs = two_species_approx_wavenumber(ωs, medium, [sp1,sp2]);
  # Approximate challis
  kTCs = multispecies_challis(ωs, medium, [sp1,sp2]);

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
  plot(xs, ys_arr, xlabel=L"a_S k", ylabel="sound speed (m/s)", labels=labs, line = styles, xlims = (-0.002,maximum(xs))
        , ylims = ( min(ys_arr[1]...,ys_arr[3]...,ys_arr[4]...)*0.995,  max(ys_arr[1]...,ys_arr[3]...,ys_arr[4]...)*1.005));
  p1 = gray_square!([xs[1],xs[m]],[y1,y2]);

  ys_arr = atten_arr;
  labs = [L"Effective $k_{*}$" L"Approximate $k_{*LS}$" L"Approximate $k_{*C}$"]
  plot(xs, ys_arr, labels=labs, xlabel=L"a_S k", ylabel="attenuation (1/m)", line=styles, xlims = (-0.002,maximum(xs))
        , ylims = ( min(ys_arr[1]...,ys_arr[3]...)*0.995,  max(ys_arr[1]...,ys_arr[3]...)*1.005));
  y1 = min(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  y2 = max(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  p2 = gray_square!([xs[1],xs[m]],[y1,y2]);

  plot(p1,p2)
  savefig("../images/compare_concrete.png")
  savefig("../images/compare_concrete.pdf")
  gui()

## Zoomed in version
  ωs = linspace(ωs[1],ωs[m],250) # k from 0 to 1
  m = length(ωs);
  kTs  = multispecies_wavenumber(ωs, medium, [sp1,sp2]);
  kTLSs = two_species_approx_wavenumber(ωs, medium, [sp1,sp2]);
  kTCs  = multispecies_challis(ωs, medium, [sp1,sp2]);

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
  p1 = gray_square!([xs[1],xs[m]],[y1,y2]);

  ys_arr = atten_arr;
  labs = [L"Effective $k_{*}$" L"Approximate $k_{*LS}$" L"Approximate $k_{*C}$"]
  plot(xs, ys_arr, xlabel=L"a_S k", ylabel="attenuation (1/m)", labels=labs
                 , border=false, line = styles, xlims = (0,maximum(xs)));
  y1 = min(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  y2 = max(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  dy = abs(y2- y1)/120.0;
  p2 = gray_square!([xs[1],xs[m]],[y1,y2]);

  plot(p1,p2)
  savefig("../images/compare_concrete_zoom.png")
  savefig("../images/compare_concrete_zoom.pdf")

## How small do the air-pockets need to be ?
  # ω = ωfactor*real(medium.c)/4.
  # r1s = (0.002:0.001:0.2)./ωfactor
  #
  # sp2 = Specie(ρ=inclusion2.ρ ,r=r2, c=inclusion2.c, volfrac = volfrac-vol1)
  #
  # # True wavenumber
  # kTs = map(r1s) do r1
  #   sp1 = Specie(ρ=inclusion1.ρ ,r=r1, c=inclusion1.c, volfrac = vol1)
  #   multispecies_wavenumber(ω, medium, [sp1,sp2])
  # end
  # # Approximate wavenumber
  # kTLSs = map(r1s) do r1
  #   sp1 = Specie(ρ=inclusion1.ρ ,r=r1, c=inclusion1.c, volfrac = vol1)
  #   two_species_approx_wavenumber(ω, medium, [sp1,sp2])
  # end
  # # Approximate wavenumber
  # kTCs = map(r1s) do r1
  #   sp1 = Specie(ρ=inclusion1.ρ ,r=r1, c=inclusion1.c, volfrac = vol1)
  #   multispecies_challis(ω, medium, [sp1,sp2])
  # end
  #
  # speed_arr = [ ω./real(kTs), ω./real(kTLSs), ω./real(kTCs), 0.*real(kTLSs) + real(medium.c)]
  # atten_arr = imag([kTs,kTLSs,kTCs])
  #
  # styles = [:solid :dashdot :dashdot :dot]
  # labs = [L"Effective $k_{*}$" L"Approximate $k_{*LS}$" L"Approximate $k_{*C}$" "Lime Stone"]
  # ys_arr = speed_arr;
  # xs = r1s.*ω./real(medium.c);
  # p1 = plot(xs, ys_arr, xlabel=L"a_S k_0", ylabel="sound speed", labels=labs
  #                , line = styles , ylims = (min(0,(minimum.(ys_arr))...), maximum(ys_arr[1])*1.2));
  #
  # ys_arr = atten_arr;
  # labs = [L"Effective $k_{*}$" L"Approximate $k_{*LS}$" L"Approximate $k_{*C}$"]
  # p2 = plot(xs, ys_arr, labels=labs, xlabel=L"a_S k_0", ylabel="attenuation"
  #                 , ylims = (max(-1,(minimum.(ys_arr))...), maximum(ys_arr[1])*1.04));
  # plot(p1,p2)
  # savefig("../images/compare_concrete_small_air.png")
  # savefig("../images/compare_concrete_small_air.pdf")

Plots.scalefontsizes(1/1.7)
