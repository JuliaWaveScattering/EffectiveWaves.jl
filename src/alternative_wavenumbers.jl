"wavenumber from Challis, R. E., et al. Ultrasound techniques for characterizing colloidal dispersions. Reports on progress in physics 68.7 (2005): 1541."
wavenumber_challis{T<:Number}(ωs::AbstractArray,medium::Medium, species::Array{Specie{T}}; kws...) = [wavenumber_challis(ω, medium, species; kws...) for ω in ωs]

# will sum Hankel orders only up to n =2 as in the paper "Ultrasound techniques for characterizing colloidal dispersions"
function wavenumber_challis{T}(ω::Number, medium::Medium{T}, species::Array{Specie{T}};
    radius_multiplier = 1.005, verbose = false, hankel_order=2)
  # rel_tol=0.004; radius_multiplier = 1.005
  # a12 = radius_multiplier(a1+a2)

  volume_fraction = sum(pi*sp.r^2.0*sp.num_density for sp in species)
  if volume_fraction >= 0.4
    warn("the volume fraction $(volume_fraction) is a bit too high, expect a relative error of approximately $(volume_fraction^3.0)")
  end

  kT2 = (ω/medium.c)^2.0
  kT2 += 4.0im*sum(sp.num_density*Zn(ω,sp,medium,m) for sp in species
  , m in -hankel_order:hankel_order)

  # second order number fraction, sum up too same hankel order
  kT2 += -(8.0/(ω/medium.c)^2)*sum(
    sp1.num_density^2*abs(m2-m1)*Zn(ω,sp1,medium,m1)*Zn(ω,sp1,medium,m2)
  for sp1 in species, m1 = -hankel_order:hankel_order, m2 = -hankel_order:hankel_order)

  return sqrt(kT2)
end


"low volfrac and low wavenumber. Fails badly for strong scatterers"
function one_species_low_wavenumber(ωs, medium, sp)
  φ = pi*sp.r^2*sp.num_density
  kS = 1.0/sp.c
  k = 1.0/medium.c
  P = 1 - kS^2*medium.ρ/(k^2*medium.ρ)
  Q = (medium.ρ - sp.ρ)/(medium.ρ + sp.ρ)
  (ωs./medium.c).*(1 - (φ/2)*(P+2*Q)-(φ^2/8)*(2*P^2-(P+2*Q)^2))
end

function wavenumber_far_field_low_volfrac{T}(ω::Complex{T}, medium::Medium{T}, species::Array{Specie{T}};
    rel_tol=0.0002, radius_multiplier = 1.005, verbose = false)
  # rel_tol=0.004; radius_multiplier = 1.005
  # a12 = radius_multiplier(a1+a2)

  volume_fraction = sum(pi*sp.r^2.0*sp.num_density for sp in species)
  if volume_fraction >= 0.4
    warn("the volume fraction $(volume_fraction) is a bit too high, expect a relative error of approximately $(volume_fraction^3.0)")
  end
  kT2 = (ω/medium.c)^2.0
  next_order = 4.0im*sum(sp.num_density*Zn(ω,sp,medium,0) for sp in species)
  hankel_order=1

  # sum more hankel orders until the relative error < rel_tol
  while abs(next_order/kT2) > rel_tol
    kT2 += next_order
    hankel_order +=1
    next_order = 4.0im*sum(sp.num_density*Zn(ω,sp,medium,m) for sp in species, m in (-hankel_order,hankel_order))
  end
  kT2 += next_order
  if verbose println("max Hankel order = $hankel_order") end

  # second order number fraction, sum up too same hankel order
  kT2 += -8.0*sum(
    sp1.num_density*sp2.num_density*abs(m2-m1)/(ω/medium.c)^2*Zn(ω,sp1,medium,m1)*Zn(ω,sp2,medium,m2)
  for sp1 in species, sp2 in species, m1 = -hankel_order:hankel_order, m2 = -hankel_order:hankel_order)

  return sqrt(kT2)
end
