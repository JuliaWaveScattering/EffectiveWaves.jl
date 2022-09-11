function reflection_coefficient(ω::T, dwave::DiscretePlaneWaveMode{T},
    source::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{2,Halfspace{T,2}};
     # medium::Acoustic{T,2}, specie::Specie{2}; θin::T = 0.0,
     scheme::Symbol = :trapezoidal) where T <: AbstractFloat

    k = ω / source.medium.c
    M = dwave.basis_order
    σ = k .* integration_scheme(dwave.x; scheme=scheme) # multiple by k to be the same as the non-dimensional version

    θin = transmission_angle(source,material)

    R = T(2)*number_density(material.microstructure.species) / (cos(θin)*k^2)*sum(
        im^T(m)*exp(im*k*dwave.x[j]*cos(θin) - im*θin*m)*dwave.amplitudes[j,m+M+1,1]*σ[j]
    for m=-M:M, j in eachindex(dwave.x))

    return R
end
