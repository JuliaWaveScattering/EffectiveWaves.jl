"wavenumber from Challis, R. E., et al. Ultrasound techniques for characterizing colloidal dispersions. Reports on progress in physics 68.7 (2005): 1541."
wavenumber_challis(ωs::AbstractArray,medium::Acoustic{T,2}, species::Array{Specie{T,2}}; kws...) where T<:AbstractFloat = [wavenumber_challis(ω, medium, species; kws...) for ω in ωs]

# will sum Hankel orders only up to n =2 as in the paper "Ultrasound techniques for characterizing colloidal dispersions"
function wavenumber_challis(ω::Number, medium::Acoustic{T,2}, species::Array{Specie{T,2}};
    radius_multiplier = 1.005, verbose = false, basis_order=2) where T<:AbstractFloat
  # tol=0.004; radius_multiplier = 1.005
  # a12 = radius_multiplier(a1+a2)

  volume_fraction = sum(sp.volume_fraction for sp in species)
  if volume_fraction >= 0.4
    @warn("the volume fraction $(round(100*volume_fraction))% is a bit too high, expect a relative error of approximately $(round(100*volume_fraction^3.0))%")
  end

  kT2 = (ω/medium.c)^2.0
  kT2 += 4.0im*sum(number_density(sp) * Zn(ω,sp,medium,m) for sp in species
  , m in -basis_order:basis_order)

  # second order number fraction, sum up too same hankel order
  kT2 += -(8.0/(ω/medium.c)^2)*sum(
    number_density(sp)^2 * abs(m2-m1)*Zn(ω,sp1,medium,m1)*Zn(ω,sp1,medium,m2)
  for sp1 in species, m1 = -basis_order:basis_order, m2 = -basis_order:basis_order)

  return sqrt(kT2)
end


"low volfrac and low wavenumber. Fails badly for strong scatterers"
function one_species_low_wavenumber(ωs, medium::Acoustic{T,2}, sp::Specie{T,2}) where T<:AbstractFloat
  φ = sp.volume_fraction
  kS = 1.0/sp.particle.medium.c
  k = 1.0/medium.c
  P = 1 - kS^2*medium.ρ/(k^2*medium.ρ)
  Q = (medium.ρ - sp.particle.medium.ρ)/(medium.ρ + sp.particle.medium.ρ)
  (ωs./medium.c).*(1 - (φ/2)*(P+2*Q)-(φ^2/8)*(2*P^2-(P+2*Q)^2))
end

function wavenumber_far_field_low_volfrac(ω::Complex{T}, medium::Acoustic{T,2}, species::Array{Specie{T,2}};
        tol=0.0002, #radius_multiplier = 1.005,
        verbose = false) where T<:AbstractFloat
  # tol=0.004; radius_multiplier = 1.005
  # a12 = radius_multiplier(a1+a2)

  volume_fraction = sum(sp.volume_fraction for sp in species)
  if volume_fraction >= 0.4
    @warn("the volume fraction $(round(100*volume_fraction))% is a bit too high, expect a relative error of approximately $(round(100*volume_fraction^3.0))%")
  end
  kT2 = (ω/medium.c)^2.0
  next_order = 4.0im*sum(number_density(sp) * Zn(ω,sp,medium,0) for sp in species)
  basis_order=1

  # sum more hankel orders until the relative error < tol
  while abs(next_order/kT2) > tol
    kT2 += next_order
    basis_order +=1
    next_order = 4.0im*sum(number_density(sp) * Zn(ω,sp,medium,m) for sp in species, m in (-basis_order,basis_order))
  end
  kT2 += next_order
  if verbose println("max Hankel order = $basis_order") end

  Zns_vec = [- t_matrix(s, medium, ω, basis_order) for s in species]

  # second order number fraction, sum up too same hankel order
  kT2 += -8.0*sum(
    number_density(species[l1]) * number_density(species[l2]) * abs(m2-m1) / (ω/medium.c)^2 * Zns_vec[l1][m1+basis_order+1,m1+basis_order+1] * Zns_vec[l2][m2+basis_order+1,m2+basis_order+1]
  for l1 = 1:length(species), l2 = 1:length(species), m1 = -basis_order:basis_order, m2 = -basis_order:basis_order)

  return sqrt(kT2)
end
