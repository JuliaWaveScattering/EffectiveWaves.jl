# This example is used to generate plots in the paper, "Reflection from a multi-species material and its transmitted effective wavenumber." Proc. R. Soc. (2018): 20170864.

# Everything related to ploting has been commented so that this example can run with requiring the Plots package. That way, this example can be checked everytime the package is modified.

using EffectiveWaves

# using LaTeXStrings
# using Plots
# height=500
 # pyplot(linewidth=3, size=(2.6*height,height), border=false)

 # Plots.scalefontsizes(1.7)


filename="concrete"
mediumname = "Lime stone"
## choose material

# concrete
  medium = LimeStone
  inclusion1 = AirDry
  inclusion2 = Brick

  ωfactor = 50.0;
  ω = ωfactor*real(medium.c)/2.0

  vols = 0.01:0.003:0.4
  vol_prop = 0.2
  r1 = 0.1/ωfactor; r2 = 1.0/ωfactor;

  # True wavenumber
  kTs = map(vols) do v
    sp1 = Specie(ρ=inclusion1.ρ ,r=r1, c=inclusion1.c, volfrac = v*vol_prop)
    sp2 = Specie(ρ=inclusion2.ρ ,r=r2, c=inclusion2.c, volfrac = v*(1.0 - vol_prop))
    wavenumber_low_volfrac(ω, medium, [sp1,sp2])
  end
  # low volfrac wavenumber
  kT2s = map(vols) do v
    sp1 = Specie(ρ=inclusion1.ρ ,r=r1, c=inclusion1.c, volfrac = v*vol_prop)
    sp2 = Specie(ρ=inclusion2.ρ ,r=r2, c=inclusion2.c, volfrac = v*(1.0 - vol_prop))
    wavenumber_very_low_volfrac(ω, medium, [sp1,sp2])
  end

  # Approximate wavenumber
  kTLSs = map(vols) do v
    sp1 = Specie(ρ=inclusion1.ρ ,r=r1, c=inclusion1.c, volfrac = v*vol_prop)
    sp2 = Specie(ρ=inclusion2.ρ ,r=r2, c=inclusion2.c, volfrac = v*(1.0 - vol_prop))
    two_species_approx_wavenumber(ω, medium, [sp1,sp2])
  end

  # Approximate challis
  kTCs = map(vols) do v
    sp1 = Specie(ρ=inclusion1.ρ ,r=r1, c=inclusion1.c, volfrac = v*vol_prop)
    sp2 = Specie(ρ=inclusion2.ρ ,r=r2, c=inclusion2.c, volfrac = v*(1.0 - vol_prop))
    wavenumber_challis(ω, medium, [sp1,sp2])
  end

  kTs_arr = [kTs,kT2s, kTLSs, kTCs];
  speed_arr = [ ω./real(ks) for ks in kTs_arr];
  push!(speed_arr, 0 .* real(kTs) .+ real(medium.c))
  atten_arr = imag.(kTs_arr)

  styles = [:solid :dash :dashdot :dashdot :dot]
  # labs = [L"k_{*}" L"k_{*0}" L"k_{*LS}" L"k_{*C}" mediumname]
  # xlabel = "Total volume fraction"
  ys_arr = speed_arr;
  xs = vols;
  m =5;
  # p1 = plot(xs, ys_arr, xlabel=xlabel, ylabel="sound speed (m/s)", labels=labs
                #  , line = styles )

  ys_arr = atten_arr;
  styles = [:solid :dash :dashdot :dashdot]
  # labs = [L"k_{*}" L"k_{*0}" L"k_{*LS}" L"k_{*C}"]
  # p2 = plot(xs, ys_arr, labels=labs, xlabel=xlabel, ylabel="attenuation (1/m)"
                    #   , line = styles )
                #  , ylims = (minimum(ys_arr[1]), maximum(ys_arr[1])*1.04));
  # plot(p1,p2)
  # gui()
#   try mkdir("media") end
#   savefig("media/compare_$(filename)_volfrac.png")
#   savefig("media/compare_$(filename)_volfrac.pdf")
#
# Plots.scalefontsizes(1/1.7)
