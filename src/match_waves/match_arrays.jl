"
Returns (LT, ER, (im*k^2*inv_w).*invV*conj(w_vec)), which connect the effective and average wave through α = LT*A + (im*k^2*inv_w).*invV*conj(w_vec).
The matching region is X[L:end].
"
function match_arrays(ω::T, wave_effs::Vector{EffectivePlaneWaveMode{T,2}}, L::Int, X::AbstractVector{T}, source::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{2,Halfspace{T,2}};
        scheme::Symbol = :trapezoidal) where T<:Number

    species = material.species

    a12k = T(2)*real(species[1].exclusion_distance * outer_radius(species[1]) * ω / source.medium.c)

    θin = transmission_angle(source,material)
    cos_in = dot(-conj(material.shape.normal), source.direction)

    θ_effs = [transmission_angle(w,material) for w in wave_effs]
    kcos_effs = [
        dot(-conj(material.shape.normal), w.wavevector)
    for w in wave_effs]

    J = length(X) - 1
    k = ω / source.medium.c
    hos = union(w.basis_order for w in wave_effs)

    if length(hos) > 1
        @warn("Expected all effective waves to have the same hankel order!")
    end
    if abs(1 - (X[2]-X[1])*J/X[J+1]) > eps(T)
        error("Expected a uniform mesh X starting at 0.")
    end

    ho = minimum(hos)
    S = length(species)

    t_vecs = get_t_matrices(source.medium, species, ω, ho)

    σ = integration_scheme(X[1:L]; scheme=scheme) # integration scheme: trapezoidal

    w_vec = T(2) .*
        [
            sum(
                exp(im*m*(θin - θ_effs[i]) + im*X[L]*(kcos_effs[i] - cos_in)) *
                number_density(species[s]) * wave_effs[i].amplitudes[m+ho+1,s]
            for m = -ho:ho, s = 1:S) / (cos_in * (kcos_effs[i] - cos_in))
        for i in eachindex(wave_effs)]

    G_arr = [
        (j > L) ? zero(Complex{T}) :
            T(2) * (-im)^T(m-1) * number_density(species[s]) * exp(im*m*θin - im*X[j] * cos_in) * σ[j] / cos_in
    for j = 1:(J+1), m = -ho:ho, s = 1:S]

    avg_wave_effs = [DiscretePlaneWaveMode(X, wave, material.shape) for wave in wave_effs]
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
                exp(im*X[J+1] * kcos_effs[p] - im*n*θ_effs[p]) *
                 wave_effs[p].amplitudes[n+ho+1] / (kcos_effs[p] + cos_in)
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
    data = [
        integration_scheme(XR[J:l+q]; scheme=scheme)
    for l = (J-q+1):J]

OffsetArray(data[end], J:(J+q))

schs =  map((J-q+1):J) do l
    sch = integration_scheme(XR[J:l+q]; scheme=scheme)
    # println("l:$l")
    # OffsetArray(sch, J:(l+q))
end

    σs = OffsetArray(
        [
            OffsetArray(integration_scheme(XR[J:l+q]; scheme=scheme), J:(l+q))
        for l = (J-q+1):J]
    , (J-q+1):J)

    Rs = [
        [
        (l+q <= J) ?
            zero(Complex{T}) :
            sum(
                (number_density(species[s])/(k^2)) * t_vecs[s][m+ho+1,m+ho+1] * im^T(n) * wave_effs[p].amplitudes[n+ho+1] *
                exp(im*XR[j] * kcos_effs[p] - im*n*θ_effs[p]) *
                (B_mat[j-l,n-m] - S_mat[j-l,n-m]) * σs[l][j]
            for j = J:(l+q), n = -ho:ho, s = 1:S)
        for l = 0:J, p in eachindex(wave_effs)]
    for m = -ho:ho]

    ER_mat = vcat((Es+Rs)...)

    return (LT_mat, ER_mat, (im*k^2*inv_w).*invV*conj(w_vec))
end
