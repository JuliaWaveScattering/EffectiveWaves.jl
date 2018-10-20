
"
Returns YT, which connects the effective and average wave through α = YT*A.
The matching region is X[XL:end].
"
function match_only_arrays(ω::T, wave_effs::Vector{EffectiveWave{T}}, XL::Int, X::AbstractVector{T}, medium::Medium{T}, species::Vector{Specie{T}}; θin::T = 0.0) where T<:Number

    XJ = length(X)
    k = ω/medium.c
    hos = union(w.hankel_order for w in wave_effs)
    if length(hos) > 1
        @warn("Expected all effective waves to have the same hankel order!")
    end
    ho = minimum(hos)
    s = length(species)

    σ = integration_scheme(X[1:XL]; scheme=:trapezoidal) # integration scheme: trapezoidal

    avg_wave_effs = [AverageWave(X, wave) for wave in wave_effs]
    vs = [
        [w.amplitudes[j,n+ho+1,1] for w in avg_wave_effs]
    for j = XL:XJ, n = -ho:ho]

    invV = inv(sum(conj(vs)[inds] * transpose(vs[inds])  for inds in eachindex(vs)))

    YT = hcat(
        [
            invV * [(j < XL) ? zero(Complex{T}) : conj(w.amplitudes[j,n+ho+1,1]) for w in avg_wave_effs]
        for j = 1:XJ, n = -ho:ho]...
    )

    return YT
end


"
returns (w_vec, q_arr), which leads to the extinction equation  sum(w_vec.*α) = im*k^2 + sum(q_arr[:].*A_avg[:]).
The index XL indicates that X[XL] is the first point in the matching region.
"
function extinc_arrays(ω::T, wave_effs::Vector{EffectiveWave{T}},
        XL::Int, X::AbstractVector{T},
        medium::Medium{T}, species::Vector{Specie{T}};
        θin::T = 0.0) where T<:Number

    k = ω/medium.c
    hos = union(w.hankel_order for w in wave_effs)
    if length(hos) > 1
        @warn("Expected all effective waves to have the same hankel order!")
    end
    ho = minimum(hos)
    S = length(species)

    XJ = length(X)
    σ = integration_scheme(X[1:XL]; scheme=:trapezoidal) # integration scheme: trapezoidal

    w_vec = (T(2)*k) .*
        [
            sum(
                exp(im*m*(θin - w.θ_eff) + im*X[XL]*(w.k_eff*cos(w.θ_eff) - k*cos(θin))/k) *
                species[l].num_density * w.amplitudes[m+ho+1,l]
            for m = -ho:ho, l = 1:S) / (cos(θin)*(w.k_eff*cos(w.θ_eff) - k*cos(θin)))
        for w in wave_effs]

    q_arr = [
        (j > XL) ? zero(Complex{T}) :
            T(2) * (-im)^T(m-1) * species[l].num_density * exp(im*m*θin - im*X[j]*cos(θin)) * σ[j] / cos(θin)
    for j = 1:XJ, m = -ho:ho, l = 1:S]

    return (w_vec, q_arr)
end
