
function two_species_approx_wavenumber(ω::Number, medium, species)
  sp1 = species[1]
  sp2 = species[2]

  if sp1.r > sp2.r
    @warn("method two_species_approx was designed for species[1] to be the smallest.")
  end
  vol = pi*sp1.num_density*sp1.r^2
  rhoeff =  medium.ρ*(medium.ρ + sp1.ρ - vol*(medium.ρ - sp1.ρ))/
  (medium.ρ + sp1.ρ + vol*(medium.ρ - sp1.ρ))

  kTS = wavenumber_low_volfrac(ω, medium, [sp1])
  mS = Medium(ρ=rhoeff, c= ω/kTS)
  wavenumber_low_volfrac(ω, mS, [sp2])
end

function two_species_approx_wavenumber(ωs::AbstractArray, medium, species)
  sp1 = species[1]
  sp2 = species[2]

  if sp1.r > sp2.r
    @warn("method two_species_approx was designed for species[1] to be the smallest.")
  end
  vol = pi*sp1.num_density*sp1.r^2
  rhoeff =  medium.ρ*(medium.ρ + sp1.ρ - vol*(medium.ρ - sp1.ρ))/
  (medium.ρ + sp1.ρ + vol*(medium.ρ - sp1.ρ))

  kTSs = wavenumber_low_volfrac(ωs, medium, [sp1])
  kTLs = [
  begin
    mS = Medium(ρ=rhoeff, c= ωs[i]/kTSs[i])
    wavenumber_low_volfrac(ωs[i], mS, [sp2])
  end
  for i in eachindex(ωs)]
end
