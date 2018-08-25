include("match_arrays.jl")

"Numerically solved the integral equation governing the average wave. Optionally can use wave_eff to approximate the wave away from the boundary."
function match_effective_waves(ω::T, medium::Medium{T}, specie::Specie{T};
        radius_multiplier::T = 1.005,
        tol::T = T(1e-5), θin::T = zero(T),
        wave_effs::Vector{EffectiveWave{T}} = [zero(EffectiveWave{T})], kws...
    ) where T<:Number

    k = real(ω/medium.c)
    a12k = T(2)*radius_multiplier*specie.r

    if maximum(abs(w.k_eff) for w in wave_effs) == zero(T)
        wave_effs = effective_waves(real(ω/medium.c), medium, [specie];
            radius_multiplier=radius_multiplier, tol=tol, θin=θin, kws...)
    end

    L, X =  X_match_waves(k, wave_effs; a12k = a12k, tol = tol)

    M = wave_effs[1].hankel_order
    J = length(collect(X))
    len = J  * (2M + 1)

    (MM_quad,b_mat) = average_wave_system(ω, X, medium, specie; tol=tol,
        radius_multiplier=radius_multiplier, hankel_order=M, θin=θin,  kws...);
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    (LT_mat, E_mat, b_eff) = match_arrays(ω, wave_effs, L, X, medium, [specie]; θin=θin);

    E_mat*LT_mat + MM_mat # double check this E_mat for where l should run

    As = MM_mat\b
    As_mat = reshape(As, (J, 2M+1, 1))

    return AverageWave(M, collect(X), As_mat)
end

"Returns (X,L), where X[L:end] is the mesh used to match with wave_effs."
function X_match_waves(k::T, wave_effs::Vector{EffectiveWave{T}};
        tol::T = T(1e-5),  a12k::T = zero(T),
        min_X::T = (-log(tol))*k/abs(cos(wave_effs[end].θ_eff)*imag(wave_effs[end].k_eff)),
        ) where T<:AbstractFloat

    #= The default options result in:
        abs(exp(im*min_X*cos(θ_effs[end])*k_effs[end]/k)) < tol
    =#
    # Based on Simpsons error = d4f * dX^5/90
    df = maximum(abs(w.k_eff * cos(w.θ_eff) / k) for w in wave_effs[1:2])
    # Based on Simpson's rule
        # dX  = (tol*90 / (df^4))^(1/5)
    # Based on trapezoidal integration
        dX  = (tol * 24 / (df^2))^(1/3)

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
