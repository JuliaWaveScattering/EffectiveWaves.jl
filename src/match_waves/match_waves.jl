include("match_arrays.jl")

"A type for matched waves."
type MatchWave{T<:AbstractFloat}
    effective_waves::Vector{EffectiveWave{T}}
    average_wave::AverageWave{T}
    match_index::Int # waves are matched between average_wave.x[match_index:end]
end


function match_effective_waves(ω::T, medium::Medium{T}, specie::Specie{T};
        radius_multiplier::T = 1.005,
        tol::T = T(1e-5), θin::T = zero(T),
        wave_effs::Vector{EffectiveWave{T}} = [zero(EffectiveWave{T})], kws...
    ) where T<:Number

    k = real(ω/medium.c)
    a12k = T(2)*radius_multiplier*specie.r

    if maximum(abs(w.k_eff) for w in wave_effs) == zero(T)
        wave_effs = effective_waves(real(ω/medium.c), medium, [specie];
            radius_multiplier=radius_multiplier,
            extinction_rescale = false,
            tol=tol, θin=θin, kws...)
    end

    L, X =  X_match_waves(k, wave_effs; a12k = a12k, tol = tol);

    # avg_wave_effs = [AverageWave(real(k), wave, X[L-1:L+1]) for wave in wave_effs]
    # ref_norm = norm(avg_wave_effs[1].amplitudes[end,:,:])/norm(wave_effs[1].amplitudes)
    # find(
    #     norm(avg_wave_effs[i].amplitudes[end,:,:])/norm(wave_effs[i].amplitudes) < ref_norm*tol
    # for i in eachindex(wave_effs))

    M = wave_effs[1].hankel_order
    J = length(collect(X))
    len = J  * (2M + 1)

    (MM_quad,b_mat) = average_wave_system(ω, X, medium, specie; tol=tol,
        radius_multiplier=radius_multiplier, hankel_order=M, θin=θin,  kws...);
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    (LT_mat, E_mat, b_eff) = match_arrays(ω, wave_effs, L, X, medium, [specie]; θin=θin);

    B = b - E_mat*b_eff
    As = (E_mat*LT_mat + MM_mat)\B
    As_mat = reshape(As, (J, 2M+1, 1))

    αs = LT_mat*As + b_eff
    # use these αs to correct the magnitude of the amplitudes of the effective waves
    for i in eachindex(wave_effs)
        wave_effs[i].amplitudes = αs[i] .* wave_effs[i].amplitudes
    end

    match_wave = MatchWave(
        wave_effs,
        AverageWave(M, collect(X), As_mat),
        L)

    return match_waves
end

"Returns (X,L), where X[L:end] is the mesh used to match with wave_effs."
function X_match_waves(k::T, wave_effs::Vector{EffectiveWave{T}};
        tol::T = T(1e-5),  a12k::T = zero(T),
        min_X::T = (-log(tol))*k/abs(cos(wave_effs[end].θ_eff)*imag(wave_effs[end].k_eff)),
        ) where T<:AbstractFloat

    #= The default min_X result in:
        abs(exp(im*min_X*cos(θ_effs[end])*k_effs[end]/k)) < tol
    =#
    # estimate a reasonable derivative. Using df2 would be too large!
    df1 = abs(wave_effs[1].k_eff * cos(wave_effs[1].θ_eff) / k)
    # take the next index, if there is one
    i =  mod(1,length(wave_effs)) + 1
    df2 = abs(wave_effs[i].k_eff * cos(wave_effs[i].θ_eff) / k)
    df = (df1+df2)/T(2)

    # Based on Simpson's rule
        # dX  = (tol*90 / (df^4))^(1/5)
    # Based on trapezoidal integration
        dX  = (tol * 24 / (df1^2))^(1/3)

    # if whole correct size a12k was given, then make dX/a12k = integer
    if a12k  != zero(T)
        n = ceil(a12k / dX)
        dX = a12k/n
    end
    max_X = min_X + 2*length(wave_effs)*dX # add points in matching region
    X = 0:dX:max_X
    L = findmin(abs.(X .- min_X))[2]

    return L, X
end
