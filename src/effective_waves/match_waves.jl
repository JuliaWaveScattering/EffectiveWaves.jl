"
Returns (L, E, (im*k^2*inv_w).*invV*w_vec), which connect the effective and average wave through α = L*A + (im*k^2*inv_w).*invV*w_vec.
The matching region is X[XL:end].
"
function match_arrays(ω::T, wave_effs::Vector{EffectiveWave{T}}, XL::Int, X::AbstractVector{T}, medium::Medium{T}, species::Vector{Specie{T}}; θin::T = 0.0) where T<:Number

    XJ = length(X)
    k = ω/medium.c
    hos = union(w.hankel_order for w in wave_effs)
    if length(hos) > 1
        warn("Expected all effective waves to have the same hankel order!")
    end
    ho = minimum(hos)
    s = length(species)

    Z = OffsetArray{Complex{Float64}}(-ho:ho, 1:s);
    for m = 0:ho, l = 1:s
        Z[m,l] = Zn(ω,species[l],medium,m)
        Z[-m,l] = Z[m,l]
    end

    σ = integration_scheme(X[1:XL]; scheme=:trapezoidal) # integration scheme: trapezoidal

    w_vec = (T(2)*k) .*
        [
            sum(
                exp(im*m*(θin - w.θ_eff) + im*X[XL]*(w.k_eff*cos(w.θ_eff) - k*cos(θin))/k) *
                species[l].num_density * Z[m,l] *  w.amplitudes[m+ho+1,l]
            for m = -ho:ho, l = 1:s) / (cos(θin)*(w.k_eff*cos(w.θ_eff) - k*cos(θin)))
        for w in wave_effs]

    q_arr = [
        (j > XL) ? zero(Complex{T}) :
            T(2) * (-im)^T(m-1) * species[l].num_density * exp(im*m*θin - im*X[j]*cos(θin)) * Z[m,l] * σ[j] / cos(θin)
    for j = 1:XJ, m = -ho:ho, l = 1:s]

    avg_wave_effs = [AverageWave(real(k), wave, X) for wave in wave_effs]
    vs = [
        [w.amplitudes[j,n+ho+1,1] for w in avg_wave_effs]
    for j = XL:XJ, n = -ho:ho]

    invV = inv(sum(conj(vs)[inds] * transpose(vs[inds])  for inds in eachindex(vs)))
    inv_w = one(T)/(transpose(w_vec)*invV*conj(w_vec))

    invVY = hcat(
        [
            invV * [(j < XL) ? zero(Complex{T}) : conj(w.amplitudes[j,n+ho+1,1]) for w in avg_wave_effs]
        for j = 1:XJ, n = -ho:ho]...
    )

    L_mat = invVY + invV * (w_vec.*inv_w) * (transpose(q_arr[:]) - transpose(w_vec)*invVY)

    S_mat = OffsetArray{Complex{Float64}}(1:XJ, -2ho:2ho);
    for j = 1:XJ, m = -2ho:2ho
        S_mat[j,m] = integrate_S(m, X[XJ] - X[j]; θin = θin)
    end

    Es = k *[
        [
            sum(
                species[l].num_density * Z[n,l] * im^T(n+1) * S_mat[j,n-m] *
                exp(im*X[XJ]*wave_effs[p].k_eff*cos(wave_effs[p].θ_eff)/k - im*n*wave_effs[p].θ_eff) *
                 wave_effs[p].amplitudes[n+ho+1] / (wave_effs[p].k_eff*cos(wave_effs[p].θ_eff) + k*cos(θin))
            for n = -ho:ho, l = 1:s)
        for j = 1:XL, p in eachindex(wave_effs)]
    for m = -ho:ho]
    E_mat = vcat(Es...)

    return (L_mat, E_mat, (im*k^2*inv_w).*invV*w_vec)
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
        warn("Expected all effective waves to have the same hankel order!")
    end
    ho = minimum(hos)
    S = length(species)

    Z = OffsetArray{Complex{Float64}}(-ho:ho, 1:S);
    for m = 0:ho, l = 1:S
        Z[m,l] = Zn(ω,species[l],medium,m)
        Z[-m,l] = Z[m,l]
    end

    XJ = length(X)
    σ = integration_scheme(X[1:XL]; scheme=:trapezoidal) # integration scheme: trapezoidal

    w_vec = (T(2)*k) .*
        [
            sum(
                exp(im*m*(θin - w.θ_eff) + im*X[XL]*(w.k_eff*cos(w.θ_eff) - k*cos(θin))/k) *
                species[l].num_density * Z[m,l] *  w.amplitudes[m+ho+1,l]
            for m = -ho:ho, l = 1:S) / (cos(θin)*(w.k_eff*cos(w.θ_eff) - k*cos(θin)))
        for w in wave_effs]

    q_arr = [
        (j > XL) ? zero(Complex{T}) :
            T(2) * (-im)^T(m-1) * species[l].num_density * exp(im*m*θin - im*X[j]*cos(θin)) * Z[m,l] * σ[j] / cos(θin)
    for j = 1:XJ, m = -ho:ho, l = 1:S]

    return (w_vec, q_arr)
end
