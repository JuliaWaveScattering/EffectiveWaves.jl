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

function kernelW3D(k::Union{T,Complex{T}}, keff::Complex{T}, pair_rs::AbstractVector{T}, gs::Matrix{V}, hks::AbstractVector{W}, basis_order::Int) where {T<:AbstractFloat, V<:AbstractVector{T}, W<:AbstractVector}

    S = size(gs,1)

    jkeffs = [
        sbesselj.(l, keff .* pair_rs)
    for l in 0:(2basis_order+1)]

    Wkers = [
        (keff .* hks[l+1] .* jkeffs[l+2] - k .* hks[l+2] .* jkeffs[l+1]) .* pair_rs.^2
    for l = 0:2basis_order]

    Wkers2 = [
        (circshift(Wkers[l],-1) - Wkers[l])[1:end-1]
    for l in eachindex(Wkers)]

    Ws = [
        sum(Wkers2[l+1] .* gs[i])
    for l = 0:2basis_order, i in CartesianIndices(gs)]

    return Ws
end

function kernelW2D(k::Union{T,Complex{T}}, keff::Complex{T}, pair_rs::AbstractVector{T}, gs::Matrix{V}, hks::AbstractVector{W}, basis_order::Int) where {T<:AbstractFloat, V<:AbstractVector{T}, W<:AbstractVector}

    S = size(gs,1)

    jkeffs = [
        besselj.(m, keff .* pair_rs)
    for m in -2basis_order:(1+2basis_order)]

    Wkers = [
        (keff .* hks[l+1] .* jkeffs[l+2] - k .* hks[l+2] .* jkeffs[l+1]) .* pair_rs
    for l = 0:4basis_order]

    Wkers2 = [
        (circshift(Wkers[l],-1) - Wkers[l])[1:end-1]
    for l in eachindex(Wkers)]

    Ws = [
        sum(Wkers2[l+1] .* gs[i])
    for l = 0:4basis_order, i in CartesianIndices(gs)]

    return Ws
end

function precalculate_pair_correlations(micro::Microstructure{3}, k::Union{T,Complex{T}} where T, basis_order::Int)

    pair_rs = micro.paircorrelations[1].r

    hks =  [shankelh1.(l, k .* pair_rs) for l in 0:(2basis_order+1)]

    # calculate average of dp between pair_rs[j] and pair_rs[j+1]
    gs = map(micro.paircorrelations) do p
        (circshift(p.dp,-1) + p.dp)[1:end-1] ./ 2
    end

    return pair_rs, hks, gs
end

function precalculate_pair_correlations(micro::Microstructure{2}, k::Union{T,Complex{T}} where T, basis_order::Int)

    pair_rs = micro.paircorrelations[1].r

    hks = [hankelh1.(m, k .* pair_rs) for m in -2basis_order:(1+2basis_order)]

    gs = map(micro.paircorrelations) do p
        # calculate segments of integrals between r_j and r_j+1
        (circshift(p.dp,-1) + p.dp)[1:end-1] ./ 2
    end

    return pair_rs, hks, gs
end
