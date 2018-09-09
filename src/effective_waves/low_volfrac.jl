wavenumber_low_volfrac(ω::Number, medium::Medium, specie::Specie; kws...) =
  wavenumber_low_volfrac(ω, medium, [specie]; kws...)

wavenumber_low_volfrac(ωs::AbstractVector{T}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number =
  [wavenumber_low_volfrac(ω, medium, species; kws...) for ω in ωs]

wavenumber_low_volfrac(ωs::AbstractVector{T}, medium::Medium{T}, specie::Specie{T}; kws...) where T<:Number = [wavenumber_low_volfrac(ω, medium, [specie]; kws...) for ω in ωs]

function wavenumber_low_volfrac(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; tol::T =1e-6,
    hankel_order::Int = maximum_hankel_order(ω, medium, species; tol=tol),
    radius_multiplier::T = 1.005, verbose::Bool = true) where T <: Number

  volume_fraction = sum(pi*sp.r^2.0*sp.num_density for sp in species)
  if volume_fraction >= 0.4 && verbose
    warn("the volume fraction $(round(100*volume_fraction))% is too high, expect a relative error of approximately $(round(100*volume_fraction^3.0))%")
  end
  num_density = sum(sp.num_density for sp in species)
  # Add incident wavenumber
  kT2 = (ω/medium.c)^2.0
  # Add far-field contribution
  kT2 += - 4.0im*num_density*far_field_pattern(ω, medium, species; hankel_order=hankel_order)(0.0)
  # Add pair-field contribution
  kT2 += - 4.0im*num_density^(2.0)*pair_field_pattern(ω, medium, species;
      radius_multiplier=radius_multiplier, hankel_order=hankel_order)(0.0)

  return (imag(sqrt(kT2)) > zero(T)) ? sqrt(kT2) : -sqrt(kT2)
end

reflection_coefficient_low_volfrac(ωs::AbstractVector{T}, medium::Medium{T}, species::Vector{Specie{T}}; kws... ) where T<:Number =
    [reflection_coefficient_low_volfrac(ω, medium, species; kws... ) for ω in ωs]

function reflection_coefficient_low_volfrac(ω::T, medium::Medium{T}, species::Vector{Specie{T}};
        θin::T = zero(T), kws... ) where T<:Number

    θ_ref = T(π) - T(2)*θin
    fo = far_field_pattern(ω, medium, species; kws...)
    dfo = diff_far_field_pattern(ω, medium, species; kws...)
    foo = pair_field_pattern(ω, medium, species; kws...)

    k = ω/medium.c
    α = k*cos(θin)
    R1 = im*fo(θ_ref)
    R2 = 2.0*fo(zero(T))/(α^2.0)
    R2 = im*foo(θ_ref) + R2*(sin(θin)*cos(θin)*dfo(θ_ref) - fo(θ_ref))

    num_density = sum(s.num_density for s in species)
    R = (R1 + num_density*R2)*num_density/(α^2.0)
    return R
end

function wavenumber_very_low_volfrac(ω::Number, medium::Medium, species::Array{Specie{T}};
    tol=1e-6, radius_multiplier = 1.005, verbose = false)  where T<:Number

  volume_fraction = sum(pi*sp.r^2.0*sp.num_density for sp in species)
  if volume_fraction >= 0.4
    warn("the volume fraction $(round(100*volume_fraction))% is too high, expect a relative error of approximately $(round(100*volume_fraction^3.0))%")
  end
  kT2 = (ω/medium.c)^2.0
  next_order = 4.0im*sum(sp.num_density*Zn(ω,sp,medium,0) for sp in species)
  hankel_order=1

  # sum more hankel orders until the relative error < tol
  while abs(next_order/kT2) > tol
    kT2 += next_order
    next_order = 4.0im*sum(sp.num_density*Zn(ω,sp,medium,m) for sp in species, m in (-hankel_order,hankel_order))
    hankel_order +=1
  end
  kT2 += next_order
  hankel_order +=1
  if verbose println("max Hankel order = $hankel_order") end

  return sqrt(kT2)
end
