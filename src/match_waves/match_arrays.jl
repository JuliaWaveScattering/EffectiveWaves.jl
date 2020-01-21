"
Returns (LT, ER, (im*k^2*inv_w).*invV*conj(w_vec)), which connect the effective and average wave through α = LT*A + (im*k^2*inv_w).*invV*conj(w_vec).
The matching region is X[L:end].
"
function match_arrays(ω::T, wave_effs::Vector{EffectivePlaneWaveMode{T}}, L::Int, X::AbstractVector{T}, medium::Acoustic{T,2}, species::Species{T,2};
        # a12k::T = 1.005*T(2)*real(specie.r*ω/medium.c),
        scheme::Symbol = :trapezoidal, θin::T = 0.0) where T<:Number

    a12k = T(2)*real(species[1].exclusion_distance * outer_radius(species[1]) * ω/medium.c)

    J = length(X) - 1
    k = ω/medium.c
    hos = union(w.basis_order for w in wave_effs)

    if length(hos) > 1
        @warn("Expected all effective waves to have the same hankel order!")
    end
    if abs(1 - (X[2]-X[1])*J/X[J+1]) > eps(T)
        error("Expected a uniform mesh X starting at 0.")
    end

    ho = minimum(hos)
    S = length(species)

    t_vecs = get_t_matrices(medium, species, ω, ho)

    σ = integration_scheme(X[1:L]; scheme=scheme) # integration scheme: trapezoidal

    w_vec = T(2) .*
        [
            sum(
                exp(im*m*(θin - w.θ_eff) + im*X[L]*(w.k_eff*cos(w.θ_eff) - cos(θin))) *
                number_density(species[s]) *  w.amplitudes[m+ho+1,s]
            for m = -ho:ho, s = 1:S) / (cos(θin)*(w.k_eff*cos(w.θ_eff) - cos(θin)))
        for w in wave_effs]

    G_arr = [
        (j > L) ? zero(Complex{T}) :
            T(2) * (-im)^T(m-1) * number_density(species[s]) * exp(im*m*θin - im*X[j]*cos(θin)) * σ[j] / cos(θin)
    for j = 1:(J+1), m = -ho:ho, s = 1:S]

    avg_wave_effs = [DiscretePlaneWaveMode(X, wave) for wave in wave_effs]
    vs = [
        [w.amplitudes[j,n+ho+1,1] for w in avg_wave_effs]
    for j = L:(J+1), n = -ho:ho]

    invV = inv(sum(conj(vs)[inds] * transpose(vs[inds])  for inds in eachindex(vs)))
    inv_w = one(T)/(transpose(w_vec)*invV*conj(w_vec))

    YTs = [
            invV * [(j < L) ? zero(Complex{T}) : conj(w.amplitudes[j,n+ho+1,1]) for w in avg_wave_effs]
        for j = 1:(J+1), n = -ho:ho]
    YT = hcat(YTs...)

    LT_mat = YT + invV * (conj(w_vec).*inv_w) * (transpose(G_arr[:]) - transpose(w_vec)*YT)

    S_mat = OffsetArray{Complex{T}}(undef, 0:J, -2ho:2ho);
    for j = 0:J, m = -2ho:2ho
        S_mat[j,m] = integrate_S(m, X[j+1]; θin = θin)
    end
    Es = [
        [
            sum(
                (number_density(species[s])/(k^2)) * t_vecs[s][m+ho+1,m+ho+1] * im^T(n+1) * S_mat[J-l,n-m] *
                exp(im*X[J+1]*wave_effs[p].k_eff*cos(wave_effs[p].θ_eff) - im*n*wave_effs[p].θ_eff) *
                 wave_effs[p].amplitudes[n+ho+1] / (wave_effs[p].k_eff*cos(wave_effs[p].θ_eff) + cos(θin))
            for n = -ho:ho, s = 1:S)
        for l = 0:J, p in eachindex(wave_effs)]
    for m = -ho:ho]

    q = min(Int(floor(a12k/(X[2]-X[1]))), J)
    B_mat = OffsetArray{Complex{T}}(undef, 0:q, -2ho:2ho);
    for j = 0:q, m = -2ho:2ho
        B_mat[j,m] = integrate_B(m, X[j+1], sqrt(abs(a12k^2 -X[j+1]^2)); θin = θin)
    end
    XR = OffsetArray((J:(J+q))*(X[2]-X[1]), J:(J+q));
    # the integration scheme changes with the domain
    σs = OffsetArray(
        [OffsetArray(integration_scheme(XR[J:l+q]; scheme=scheme), J:(l+q)) for l = (J-q+1):J]
    , (J-q+1):J)

    Rs = [
        [
        (l+q <= J) ?
            zero(Complex{T}) :
            sum(
                (number_density(species[s])/(k^2)) * t_vecs[s][m+ho+1,m+ho+1] * im^T(n) * wave_effs[p].amplitudes[n+ho+1] *
                exp(im*XR[j]*wave_effs[p].k_eff*cos(wave_effs[p].θ_eff) - im*n*wave_effs[p].θ_eff) *
                (B_mat[j-l,n-m] - S_mat[j-l,n-m]) * σs[l][j]
            for j = J:(l+q), n = -ho:ho, s = 1:S)
        for l = 0:J, p in eachindex(wave_effs)]
    for m = -ho:ho]

    ER_mat = vcat((Es+Rs)...)

    return (LT_mat, ER_mat, (im*k^2*inv_w).*invV*conj(w_vec))
end
