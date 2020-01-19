
function two_species_approx_wavenumber(ω::Number, medium, species)
  sp1 = species[1]
  sp2 = species[2]

  if outer_radius(sp1) > outer_radius(sp2)
    @warn("method two_species_approx was designed for species[1] to be the smallest.")
  end
  vol = volume_fraction(sp1)
  rhoeff =  medium.ρ*(medium.ρ + sp1.particle.medium.ρ - vol*(medium.ρ - sp1.particle.medium.ρ)) /
  (medium.ρ + sp1.particle.medium.ρ + vol*(medium.ρ - sp1.particle.medium.ρ))

  kTS = wavenumber_low_volfrac(ω, medium, [sp1])
  mS = Acoustic(2; ρ=rhoeff, c= ω/kTS)
  wavenumber_low_volfrac(ω, mS, [sp2])
end

function two_species_approx_wavenumber(ωs::AbstractArray, medium, species)
  sp1 = species[1]
  sp2 = species[2]

  if outer_radius(sp1) > outer_radius(sp2)
    @warn("method two_species_approx was designed for species[1] to be the smallest.")
  end
  vol = volume_fraction(sp1)
  rhoeff =  medium.ρ*(medium.ρ + sp1.particle.medium.ρ - vol*(medium.ρ - sp1.particle.medium.ρ))/
  (medium.ρ + sp1.particle.medium.ρ + vol*(medium.ρ - sp1.particle.medium.ρ))

  kTSs = wavenumber_low_volfrac(ωs, medium, [sp1])
  kTLs = [
  begin
    mS = Acoustic(2; ρ=rhoeff, c= ωs[i]/kTSs[i])
    wavenumber_low_volfrac(ωs[i], mS, [sp2])
  end
  for i in eachindex(ωs)]
end
