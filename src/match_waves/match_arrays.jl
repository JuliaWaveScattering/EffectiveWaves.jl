"
Returns (LT, E, (im*k^2*inv_w).*invV*conj(w_vec)), which connect the effective and average wave through α = LT*A + (im*k^2*inv_w).*invV*conj(w_vec).
The matching region is X[L:end].
"
function match_arrays(ω::T, wave_effs::Vector{EffectiveWave{T}}, L::Int, X::AbstractVector{T},
        medium::Medium{T}, species::Vector{Specie{T}};
        a12k::T = 1.005*T(2)*real(specie.r*ω/medium.c),
        scheme::Symbol = :trapezoidal, θin::T = 0.0) where T<:Number

    J = length(X) - 1
    k = ω/medium.c
    hos = union(w.hankel_order for w in wave_effs)

    if length(hos) > 1
        warn("Expected all effective waves to have the same hankel order!")
    end
    if abs(1 - (X[2]-X[1])*J/X[J+1]) > eps(T)
        error("Expected a uniform mesh X starting at 0.")
    end

    ho = minimum(hos)
    S = length(species)

    Z = OffsetArray{Complex{Float64}}(-ho:ho, 1:S);
    for m = 0:ho, s = 1:S
        Z[m,s] = Zn(ω,species[s],medium,m)
        Z[-m,s] = Z[m,s]
    end

    σ = integration_scheme(X[1:L]; scheme=scheme) # integration scheme: trapezoidal

    w_vec = T(2) .*
        [
            sum(
                exp(im*m*(θin - w.θ_eff) + im*X[L]*(w.k_eff*cos(w.θ_eff) - cos(θin))) *
                species[s].num_density * Z[m,s] *  w.amplitudes[m+ho+1,s]
            for m = -ho:ho, s = 1:S) / (cos(θin)*(w.k_eff*cos(w.θ_eff) - cos(θin)))
        for w in wave_effs]

    G_arr = [
        (j > L) ? zero(Complex{T}) :
            T(2) * (-im)^T(m-1) * species[s].num_density * exp(im*m*θin - im*X[j]*cos(θin)) * Z[m,s] * σ[j] / cos(θin)
    for j = 1:(J+1), m = -ho:ho, s = 1:S]

    avg_wave_effs = [AverageWave(wave, X) for wave in wave_effs]
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

    S_mat = OffsetArray{Complex{T}}(0:J, -2ho:2ho);
    for j = 0:J, m = -2ho:2ho
        S_mat[j,m] = integrate_S(m, X[j+1]; θin = θin)
    end
    Es = [
        [
            sum(
                species[s].num_density * Z[n,s] * im^T(n+1) * S_mat[J-l,n-m] *
                exp(im*X[J+1]*wave_effs[p].k_eff*cos(wave_effs[p].θ_eff) - im*n*wave_effs[p].θ_eff) *
                 wave_effs[p].amplitudes[n+ho+1] / (wave_effs[p].k_eff*cos(wave_effs[p].θ_eff) + cos(θin))
            for n = -ho:ho, s = 1:S)
        for l = 0:J, p in eachindex(wave_effs)]
    for m = -ho:ho]

    # q = min(Int(floor(a12k/(X[2]-X[1]))), J)
    # B_mat = OffsetArray{Complex{T}}(0:q, -2ho:2ho);
    # for j = 0:q, m = -2ho:2ho
    #     B_mat[j,m] = integrate_B(m, X[j+1], sqrt(abs(a12k^2 -X[j]^2)); θin = θin)
    # end
    # XR = OffsetArray((J:(J+q))*(X[2]-X[1]), J:(J+q));
    # Rs = [
    #     [
    #         (l+q <= J) ?
    #             zero(Complex{T}) :
    #             begin
    #                 # the integration scheme may change with the domain
    #                 σ = OffsetArray(integration_scheme(XR[J:l+q]; scheme=scheme), J:(l+q));
    #                 sum(
    #                     species[s].num_density * Z[n,s] * im^T(n) * wave_effs[p].amplitudes[n+ho+1] *
    #                     exp(im*XR[j]*wave_effs[p].k_eff*cos(wave_effs[p].θ_eff) - im*n*wave_effs[p].θ_eff) *
    #                     (B_mat[j-l,n-m] - S_mat[j-l,n-m]) * σ[j]
    #                 for j = J:(l+q) n = -ho:ho, s = 1:S)
    #             end
    #     for l = 0:J, p in eachindex(wave_effs)]
    # for m = -ho:ho]

    E_mat = vcat(Es...)

    return (LT_mat, E_mat, (im*k^2*inv_w).*invV*conj(w_vec))
end


"
Returns YT, which connects the effective and average wave through α = YT*A.
The matching region is X[XL:end].
"
function match_only_arrays(ω::T, wave_effs::Vector{EffectiveWave{T}}, XL::Int, X::AbstractVector{T}, medium::Medium{T}, species::Vector{Specie{T}}; θin::T = 0.0) where T<:Number

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

    avg_wave_effs = [AverageWave(real(k), wave, X) for wave in wave_effs]
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
