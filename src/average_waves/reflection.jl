function reflection_coefficient(ω::T, amps::AverageWave{T}, medium::Medium, specie::Specie;
        θin::T = 0.0, scheme::Symbol = :trapezoidal) where T <: AbstractFloat

    k = ω/medium.c
    M = amps.hankel_order
    σ = k .* integration_scheme(amps.x; scheme=scheme) # multiple by k to be the same as the non-dimensional version

    R = T(2)*specie.num_density/(cos(θin)*k^2)*sum(
        im^T(m)*exp(im*k*amps.x[j]*cos(θin) - im*θin*m)*amps.amplitudes[j,m+M+1,1]*σ[j]
    for m=-M:M, j in eachindex(amps.x))

    return R
end
