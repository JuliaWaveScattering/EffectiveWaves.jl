
d(x,m) = diffbesselj(m,x)*diffhankelh1(m,x) + (1.0 - (m/x)^2)*besselj(m,x)*hankelh1(m,x)

function one_species_low_wavenumber(ωs, medium, sp)
  φ = pi*sp.r^2*sp.num_density
  kS = 1.0/sp.c
  k = 1.0/medium.c
  P = 1 - kS^2*medium.ρ/(k^2*medium.ρ)
  Q = (medium.ρ - sp.ρ)/(medium.ρ + sp.ρ)
  (ωs./medium.c).*(1 - (φ/2)*(P+2*Q)-(φ^2/8)*(2*P^2-(P+2*Q)^2))
end

wavenumber_low_volfrac(ω::Number, medium::Medium, specie::Specie; kws...) =
  wavenumber_low_volfrac(ω, medium, [specie]; kws...)

wavenumber_low_volfrac{T<:Number}(ωs::AbstractArray,medium::Medium, species::Array{Specie{T}}; kws...) =
  [wavenumber_low_volfrac(ω, medium, species; kws...) for ω in ωs]

wavenumber_low_volfrac(ωs::AbstractArray,medium::Medium, specie::Specie; kws...) =
  [wavenumber_low_volfrac(ω, medium, [specie]; kws...) for ω in ωs]

# test{T<:Number}(x,species::Array{Specie{T}}) = println(species)
# test2(x,species::Array{Specie{T}}) where {T<:Number} = prinlnt(species)

function wavenumber_low_volfrac{T<:Number}(ω::Number, medium::Medium, species::Array{Specie{T}};
    rel_tol=2e-5, radius_multiplier = 1.005, verbose = false)
  # rel_tol=0.004; radius_multiplier = 1.005
  # a12 = radius_multiplier(a1+a2)

  volume_fraction = sum(pi*sp.r^2.0*sp.num_density for sp in species)
  if volume_fraction >= 0.4
    warn("the volume fraction $(volume_fraction) is too high, expect a relative error of approximately $(volume_fraction^3.0)")
  end
  kT2 = (ω/medium.c)^2.0
  next_order = 4.0im*sum(sp.num_density*Zn(ω,sp,medium,0) for sp in species)
  hankel_order=1

  # sum more hankel orders until the relative error < rel_tol
  while abs(next_order/kT2) > rel_tol
    kT2 += next_order
    next_order = 4.0im*sum(sp.num_density*Zn(ω,sp,medium,m) for sp in species, m in (-hankel_order,hankel_order))
    hankel_order +=1
  end
  kT2 += next_order
  hankel_order +=1
  if verbose println("max Hankel order = $hankel_order") end

  # second order number fraction
  second_order(m1_arr, m2_arr) = 4.0im*pi*sum(
    begin
      a12 = radius_multiplier*(sp1.r + sp2.r)
      sp1.num_density*sp2.num_density*a12^2.0*d(ω/medium.c*a12,m2-m1)*
      Zn(ω,sp1,medium,m1)*Zn(ω,sp2,medium,m2)
    end
  for sp1 in species, sp2 in species, m1 = m1_arr, m2 = m2_arr)

  next_order += second_order(-hankel_order:hankel_order,-hankel_order:hankel_order)
  hankel_order+=1

  while abs(next_order/kT2) > rel_tol
    kT2 += next_order
    next_order = second_order(-hankel_order:hankel_order, (-hankel_order,hankel_order))
    next_order += second_order((-hankel_order,hankel_order), -(hankel_order-1):(hankel_order-1))
    hankel_order +=1
  end

  return sqrt(kT2)
end

function wavenumber_very_low_volfrac{T<:Number}(ω::Number, medium::Medium, species::Array{Specie{T}};
    rel_tol=0.0002, radius_multiplier = 1.005, verbose = false)

  volume_fraction = sum(pi*sp.r^2.0*sp.num_density for sp in species)
  if volume_fraction >= 0.4
    warn("the volume fraction $(volume_fraction) is too high, expect a relative error of approximately $(volume_fraction^3.0)")
  end
  kT2 = (ω/medium.c)^2.0
  next_order = 4.0im*sum(sp.num_density*Zn(ω,sp,medium,0) for sp in species)
  hankel_order=1

  # sum more hankel orders until the relative error < rel_tol
  while abs(next_order/kT2) > rel_tol
    kT2 += next_order
    next_order = 4.0im*sum(sp.num_density*Zn(ω,sp,medium,m) for sp in species, m in (-hankel_order,hankel_order))
    hankel_order +=1
  end
  kT2 += next_order
  hankel_order +=1
  if verbose println("max Hankel order = $hankel_order") end

  return sqrt(kT2)
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
