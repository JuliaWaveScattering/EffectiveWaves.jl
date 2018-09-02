"A type for matched waves."
type MatchWave{T<:AbstractFloat}
    effective_waves::Vector{EffectiveWave{T}}
    average_wave::AverageWave{T}
    x_match::Vector{T} # waves are matched between average_wave.x[match_index:end]
end

function MatchWave(ω::T, medium::Medium{T}, specie::Specie{T};
        radius_multiplier::T = 1.005,
        tol::T = T(1e-5), θin::T = zero(T),
        hankel_order::Int = 2,
        max_size::Int = 500,
        wave_effs::Vector{EffectiveWave{T}} = [zero(EffectiveWave{T})], kws...
    ) where T<:Number

    k = real(ω/medium.c)

    if maximum(abs(w.k_eff) for w in wave_effs) == zero(T)
        wave_effs = effective_waves(k, medium, [specie];
            radius_multiplier=radius_multiplier,
            extinction_rescale = false, #hankel_order=hankel_order,
            tol=1e-8, θin=θin,
            hankel_order=hankel_order,
            kws...)
    elseif wave_effs[1].hankel_order != hankel_order
        error("wave_effs given have a different hankel order than the option hankel_order=$(hankel_order)")
    end
    # use non-dimensional effective waves
    wave_non_effs = deepcopy(wave_effs)
    for w in wave_non_effs
       w.k_eff = w.k_eff/k
    end

    a12k = T(2)*radius_multiplier*specie.r*k
    # using non-dimensional wave_non_effs and a12k results in non-dimensional mesh X
    L, X =  x_mesh_match(wave_non_effs; a12 = a12k, tol = tol, max_size=max_size);

    avg_wave_effs = [AverageWave(X[L:L+1], w) for w in wave_non_effs]
    for i in eachindex(wave_non_effs)
        wave_non_effs[i].amplitudes = wave_non_effs[i].amplitudes / norm(avg_wave_effs[i].amplitudes[1,:,1])
        wave_effs[i].amplitudes = wave_non_effs[i].amplitudes
    end

    J = length(collect(X)) - 1
    len = (J + 1)  * (2hankel_order + 1)

    (MM_quad,b_mat) = average_wave_system(ω, X, medium, specie; tol=tol,
        radius_multiplier=radius_multiplier, hankel_order=hankel_order, θin=θin,  kws...);
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    (LT_mat, ER_mat, b_eff) = match_arrays(ω, wave_non_effs, L, X, medium, [specie]; θin=θin, a12k=a12k);

    B = b - ER_mat*b_eff
    As = (ER_mat*LT_mat + MM_mat)\B
    As_mat = reshape(As, (J+1, 2hankel_order+1, 1))

    αs = LT_mat*As + b_eff
    # use these αs to correct the magnitude of the amplitudes of the effective waves
    for i in eachindex(wave_effs)
        wave_effs[i].amplitudes = αs[i] .* wave_effs[i].amplitudes
    end

    return MatchWave(wave_effs, AverageWave(hankel_order, collect(X)./k, As_mat), collect(X[L:end])./k)
end

"Returns (x,L), where x[L:end] is the mesh used to match with wave_effs."
function x_mesh_match(wave_effs::Vector{EffectiveWave{T}}; tol::T = 1e-7, kws... ) where T<:AbstractFloat
    # below wave_effs[end] establishes how long X should be, while wave_effs[1] estalishes how fine the mesh should be.
    x = (length(wave_effs) == 1 || abs(cos(wave_effs[end].θ_eff)*imag(wave_effs[end].k_eff)) < tol*T(10) ) ?
        x_mesh(wave_effs[end], wave_effs[1]; max_x = T(2pi)/abs(wave_effs[1].k_eff), tol=tol, kws...) :
        x_mesh(wave_effs[end], wave_effs[1]; tol=tol, kws...)

    x_match = x[end]
    x_max = (length(x) < length(wave_effs)*T(1.5)) ?
         x[end] + T(1.5)*(x[2] - x[1])*length(wave_effs) :
         T(2) * x[end]

    x = 0.0:(x[2] - x[1]):x_max # choose twice the match length to heavily penalise overfitting wave_effs[end]

    L = findmin(abs.(x .- x_match))[2]

    return L, x
end
