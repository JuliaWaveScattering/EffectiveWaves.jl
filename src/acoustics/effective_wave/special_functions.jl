function kernelN2D(n::Int,x::Union{T,Complex{T}},y::Union{T,Complex{T}}) where T<:AbstractFloat
        h = hankelh1(n,x); dh = diffhankelh1(n,x)
        j = besselj(n,y);  dj = diffbesselj(n,y)

    return x * dh * j - y * h * dj
end

function kernelN3D(n::Int,x::Union{T,Complex{T}},y::Union{T,Complex{T}}) where T<:AbstractFloat
    h = shankelh1(n,x); dh = diffshankelh1(n,x)
    j = sbesselj(n,y);  dj = diffsbesselj(n,y)

    return x * dh * j - y * h * dj
end

function kernelN3D(k::Union{T,Complex{T}}, keff::Complex{T}, as::Matrix{T}, pair_rs::AbstractVector{T}, gs::Matrix{T}, hks::AbstractVector{T}, basis_order::Int) where T<:AbstractFloat

    jkeffs = [
        sbesselj.(l, keff .* pair_rs)
    for l in 0:(2basis_order+1)]

    Wkers = [
        keff .* hks[l] .* jkeffs[l+1] - k .* hks[l+1] .* jkeffs[l]
    for l = 0:2basis_order]

    Ws = [
        sum(
            (circshift(Wker[l+1],-1) - Wker[l+1])[1:end-1] .* gs[s1,s2]
        )
    for l = 0:2basis_order, s1 = 1:S, s2 = 1:S]

    Ns = [
        kernelN3D(l,k*as[s1,s2],keff*as[s1,s2])
    for l = 0:2ho, s1 = 1:S, s2 = 1:S]

    Ns = Ns + Ws;

    return Ns
end

function precalculate_pair_correlations(micro::Microstructure, basis_order)
    pair_rs = micro.paircorrelations.r
    hks = [shankelh1.(l, k .* pair_rs) for l in 0:(2basis_order+1)]
    σs = integration_scheme(pair_rs)

    gs = map(micro.paircorrelations) do p
        g = p.dp .* σs .* pair_rs .^2

        # calculate segments of integrals between r_j and r_j+1
        (circshift(g,-1) + g)[1:end-1] ./ T(2)
    end

    return pair_rs, hks, gs
end
