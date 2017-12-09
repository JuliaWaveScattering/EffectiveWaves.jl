
function two_species_approx_wavenumber(ω::Number, medium, species)
  sp1 = species[1]
  sp2 = species[2]

  if sp1.r > sp2.r
    warn("method two_species_approx was designed for species[1] to be the smallest.")
  end
  vol = pi*sp1.num_density*sp1.r^2
  rhoeff =  medium.ρ*(medium.ρ + sp1.ρ - vol*(medium.ρ - sp1.ρ))/
  (medium.ρ + sp1.ρ + vol*(medium.ρ - sp1.ρ))

  kTS = multispecies_wavenumber(ω, medium, [sp1])
  mS = Medium(ρ=rhoeff, c= ω/kTS)
  multispecies_wavenumber(ω, mS, [sp2])
end

function two_species_approx_wavenumber(ωs::AbstractArray, medium, species)
  sp1 = species[1]
  sp2 = species[2]

  if sp1.r > sp2.r
    warn("method two_species_approx was designed for species[1] to be the smallest.")
  end
  vol = pi*sp1.num_density*sp1.r^2
  rhoeff =  medium.ρ*(medium.ρ + sp1.ρ - vol*(medium.ρ - sp1.ρ))/
  (medium.ρ + sp1.ρ + vol*(medium.ρ - sp1.ρ))

  kTSs = multispecies_wavenumber(ωs, medium, [sp1])
  kTLs = [
  begin
    mS = Medium(ρ=rhoeff, c= ωs[i]/kTSs[i])
    multispecies_wavenumber(ωs[i], mS, [sp2])
  end
  for i in eachindex(ωs)]
end
