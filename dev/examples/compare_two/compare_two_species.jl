using EffectiveWaves
using LaTeXStrings
using Plots
height=500
# unicodeplots()
 pyplot(linewidth=3, size=(2.6*height,height))

 Plots.scalefontsizes(1.7)

## choose material

  # concrete
  medium = LimeStone
  inclusion1 = AirDry
  inclusion2 = Brick
  ωs = LinRange(real(medium.c/10000),real(medium.c),100) # k from 0 to 1
  # ωs = LinRange(real(medium.c/10000),real(medium.c/1000),100) # k from 0 to 1

  volfrac = 0.16
  r1 = 0.1; vol1 = 0.06
  r2 = 1.0; vol2 = volfrac - vol1

  sp1 = Specie(ρ=inclusion1.ρ ,r=r1, c=inclusion1.c, volfrac = vol1)
  sp2 = Specie(ρ=inclusion2.ρ ,r=r2, c=inclusion2.c, volfrac = volfrac-vol1)

  # True wavenumber
  kTs  = wavenumber_low_volumefraction(ωs, medium, [sp1,sp2]);
  # Approximate wavenumber
  kTLSs = two_species_approx_wavenumber(ωs, medium, [sp1,sp2]);

  speed_arr = [ ωs./real(kTs), ωs./real(kTLSs), 0 .*ωs + real(medium.c)]
  atten_arr = imag([kTs,kTLSs])

  styles = reshape([:solid,:solid,:dot],1,3)
  labs = reshape([L"Effective $k_{*}$",L"Approximate $k_{*LS}$", "Lime Stone"] ,1, 3);
  ys_arr = speed_arr;
  xs = r1.*(ωs./real(medium.c));
  m =5;
  y1 = min(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  y2 = max(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  plot(xs, ys_arr, xlabel=L"a_S k", ylabel="sound speed", labels=labs, line = styles);
  p1 = gray_square([xs[1],xs[m]],[y1,y2], plot!)

  ys_arr = atten_arr;
  labs = reshape([L"Effective $k_{*}$",L"Approximate $k_{*LS}$"] ,1, 2);
  plot(xs, ys_arr, labels=labs, xlabel=L"a_S k", ylabel="attenuation");
  y1 = min(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  y2 = max(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  p2 = gray_square([xs[1],xs[m]],[y1,y2], plot!)

  plot(p1,p2)
  savefig("media/compare_concrete.png")
  savefig("media/compare_concrete.pdf")
  gui()

## Zoomed in version
  ωs = LinRange(real(medium.c/10000),ωs[m],100) # k from 0 to 1
  m = length(ωs);
  kTs  = wavenumber_low_volumefraction(ωs, medium, [sp1,sp2]);
  kTLSs = two_species_approx_wavenumber(ωs, medium, [sp1,sp2]);

  speed_arr = [ ωs./real(kTs), ωs./real(kTLSs)]
  atten_arr = imag([kTs,kTLSs])

  labs = reshape([L"Effective $k_{*}$",L"Approximate $k_{*LS}$"] ,1, 2);
  ys_arr = speed_arr;
  xs = r1.*(ωs./real(medium.c));
  plot(xs, ys_arr, xlabel=L"a_S k", ylabel="sound speed", labels=labs
                 , border = false);
  y1 = min(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  y2 = max(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  p1 = gray_square([xs[1],xs[m]],[y1,y2], plot!)

  ys_arr = atten_arr;
  labs = reshape([L"Effective $k_{*}$",L"Approximate $k_{*LS}$"] ,1, 2);
  plot(xs, ys_arr, xlabel=L"a_S k", ylabel="attenuation", labels=labs
                 , border=false);
  y1 = min(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  y2 = max(ys_arr[1][1:m]..., ys_arr[2][1:m]...);
  dy = abs(y2- y1)/120.0;
  p2 = gray_square([xs[1],xs[m]],[y1,y2], plot!)

  plot(p1,p2)
  mkdir("media")
  savefig("media/compare_concrete_zoom.png")
  savefig("media/compare_concrete_zoom.pdf")

Plots.scalefontsizes(1/1.7)
