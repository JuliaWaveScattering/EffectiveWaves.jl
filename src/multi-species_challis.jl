try include("scattering.jl") end

multispecies_challis{T<:Number}(ωs::AbstractArray,medium::Medium, species::Array{Specie{T}}; kws...) = [multispecies_challis(ω, medium, species; kws...) for ω in ωs]


# will sum Hankel orders only up to n =2 as in the paper "Ultrasound techniques for characterizing colloidal dispersions"
function multispecies_challis{T}(ω::Number, medium::Medium{T}, species::Array{Specie{T}};
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
