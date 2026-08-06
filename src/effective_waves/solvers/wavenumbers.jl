"""
    refine_wavenumber(disp::Function, k::Complex{T}, kscale::T; tol, max_iters, max_distance)

Refine the estimate `k` of a root of the complex analytic function `disp`, such as the function returned by [`dispersion_complex`](@ref), using Newton's method with a central finite difference for the derivative.

Returns `(k_refined, converged::Bool)`. As `disp` is analytic, `disp(k_eff) = 0` is two real equations in the two real unknowns `real(k_eff)` and `imag(k_eff)`, so Newton's method converges quadratically from a nearby estimate. This is both faster and more accurate than minimising `abs(disp(k_eff))` with a derivative free method.

## Keywords
- `tol::T`: the iteration stops once a Newton step is smaller than `tol * max(abs(k), kscale)` (default: `1e-6`)
- `max_iters::Int`: the maximum number of Newton steps (default: `20`)
- `max_distance::T`: `converged` is set to `false` if an iterate moves further than this from the initial estimate, which indicates that the method has left the root being tracked (default: `kscale`)
"""
function refine_wavenumber(disp::Function, k::Complex{T}, kscale::T;
        tol::T = T(1e-6),
        max_iters::Int = 20,
        max_distance::T = kscale
    ) where T<:AbstractFloat

    k_start = k

    # a central difference is most accurate for a step of the order eps^(1/3)
    h = eps(T)^(one(T)/3) * max(abs(k_start), kscale)

    for _ in 1:max_iters
        f = disp(k)
        isfinite(f) || return (k, false)
        iszero(f) && return (k, true)

        df = (disp(k + h) - disp(k - h)) / (2*h)
        (isfinite(df) && !iszero(df)) || return (k, false)

        dk = - f / df
        isfinite(dk) || return (k, false)

        # the estimate is expected to be close to the root, so a large step means that Newton's method has been thrown off the wavenumber being tracked
        (abs(k + dk - k_start) < max_distance) || return (k, false)

        k = k + dk
        (abs(dk) < tol * max(abs(k), kscale)) && return (k, true)
    end

    return (k, false)
end

"""
    wavenumbers(ωs::AbstractVector{T}, micro::Microstructure{Dim}; symmetry, tol, branch_number, optimoptions, kws...)

Compute effective wavenumbers over a range of frequencies using adaptive interpolation and refinement.

## Arguments
- `ωs::AbstractVector{T}`: the angular frequencies, which need to be sorted in increasing order
- `tol::T`: Tolerance for numerical optimization (default: `1e-6`)
- `branch_number::Int`: Which branch of wavenumbers to return (default: `1`)
- `num_wavenumbers::Int`: How many wavenumbers to search for when using exact root-finding (default: `3`)
- `basis_order::Int`: Use the same truncation order of the basis for every frequency. When not given, the order increases with frequency in the same way as [`dispersion_complex`](@ref), but never goes below 2.
- `basis_orders::AbstractVector{Int}`: The truncation order used for each frequency in `ωs`, which overrides `basis_order`.

## Returns
- `keffs::Vector{Complex{T}}`: Complex effective wavenumbers with positive imaginary part

## Method
- The first effective wavenumber is computed at ω₁ using exact root-finding
- The second wavenumber is computed at ω₂ using exact root-finding
- `keffs[3]` is predicted using linear interpolation of `keffs[1]` and `keffs[2]`
- For i > 3, `keffs[i]` is predicted using quadratic (Lagrange) interpolation of the three previous nodes
- Each prediction is refined with [`refine_wavenumber`](@ref), falling back to a derivative free search when Newton's method does not converge

## Notes
- The method assumes small step sizes in wavenumber space (Δka < 0.001)
- The particles need to be small compared to the starting wavelength (ka < 0.1)
- Physical solutions must have Im(k_eff) > 0; negative imaginary parts are negated
- NOTE: PlanarAzimuthalSymmetry() does not include all possible wavenumbers
"""
function wavenumbers(ωs::AbstractVector{T}, micro::Microstructure{Dim};
        symmetry::AbstractSymmetry{Dim} = PlanarAzimuthalSymmetry{Dim}(),
        tol::T = 1e-6, branch_number::Int = 1,
        basis_order::Union{Int,Nothing} = nothing,
        basis_orders::AbstractVector{Int} = isnothing(basis_order) ?
            [max(2,
                3 * Int(round(maximum(outer_radius.(micro.species)) * ω / abs(micro.medium.c))) + 1
            ) for ω in ωs] :
            basis_order * ones(Int, length(ωs)),
        num_wavenumbers::Int = 3,
        bisection_method::Bool = false,
        optimoptions = Optim.Options(
            iterations = max(50, Int(round(-log(tol))) * 20),
            g_abstol = tol^T(2), x_abstol=tol^T(2)),
        kws...) where {T<:AbstractFloat,Dim}

    # check the given parameters
        isempty(ωs) && return Complex{T}[]

        if !issorted(ωs)
            throw(ArgumentError("the frequencies ωs need to be sorted in increasing order, as every wavenumber is tracked from the previous frequency."))
        end

        if branch_number < 1
            throw(ArgumentError("branch_number = $(branch_number) needs to be a positive integer."))
        end

        if length(basis_orders) != length(ωs)
            throw(ArgumentError("basis_orders has $(length(basis_orders)) elements, but needs one element for each of the $(length(ωs)) frequencies given in ωs."))
        end

    a = mean(micro.species .|> outer_radius)
    kas = a .* ωs ./ real(micro.medium.c)

    if length(kas) > 1
        dka = abs.((kas - circshift(kas,1))[2:end]) |> minimum

        if dka > 0.001
            @warn "The method is designed to take very small step increases in the wavenumber ka. The minimum step increase in ka is $(dka), but should be smaller than 0.001."
        end
    end

    if abs(kas[1]) > 0.1
        @warn "The method is designed to start with particles that are small compared to the wavelength. Starting with ka = $(kas[1]) might be too large. "
    end

    # ceff = effective_medium(micro).c
    disp = dispersion_complex(ωs[1], micro, symmetry; basis_order = basis_orders[1], kws...)

    k0s = wavenumbers(ωs[1], micro;
        symmetry = symmetry,
        basis_order = basis_orders[1],
        num_wavenumbers = num_wavenumbers,
        bisection_method = bisection_method,
        tol = tol, kws...)
    # k0s = wavenumbers_bisection_robust(ωs[1], micro; symmetry = symmetry, num_wavenumbers = 10, tol = tol, kws...)

    if isempty(k0s)
        error("No effective wavenumbers were found for the first frequency ω = $(ωs[1]). Without a starting wavenumber the other frequencies can not be tracked.")
    end

    if branch_number > length(k0s)
        @warn "Branch number $(branch_number) is larger than the number of wavenumbers found $(length(k0s)). Returning the last wavenumber found."
        branch_number = length(k0s)
    end

    errors = Vector{T}(undef, length(ωs))

    keffs = Vector{Complex{T}}(undef, length(ωs))
    keffs[1] = k0s[branch_number]
    errors[1] = disp(keffs[1]) |> abs

    interpolation_orders = zeros(Int, length(ωs))
    interpolation_orders[1] = 0

    # refine the prediction keffs[i] so that it satisfies the dispersion equation. Here kscale is the wavenumber of the background medium, which sets the scale for both the search and the convergence criteria.
    function refine!(i::Int, disp::Function, kscale::T)

        k, converged = refine_wavenumber(disp, keffs[i], kscale;
            tol = tol, max_distance = kscale)

        if !converged
            # Newton's method was thrown off, so fall back to a derivative free search. The simplex is scaled to the background wavenumber, rather than to the size of keffs[i], as otherwise a weakly attenuating wave is given a simplex which is far larger than imag(keffs[i]).
            dk = sqrt(tol) * kscale
            inner_optimizer = NelderMead(
                initial_simplex = MySimplexer(dk,dk),
                parameters = NelderMeadparameters()
            )

            f_vec(x_vec) = abs(disp(x_vec[1] + x_vec[2]*im))
            res = optimize(f_vec, [keffs[i] |> real, keffs[i] |> imag], inner_optimizer, optimoptions)

            k = res.minimizer[1] + res.minimizer[2]*im
        end

        keffs[i] = k
        errors[i] = disp(k) |> abs

        return nothing
    end

    function predict_next!(interpolation_order::Int, i::Int; kws...)

        disp = dispersion_complex(ωs[i], micro, symmetry; kws...)
        kscale = ωs[i] / real(micro.medium.c)

        if interpolation_order == 0
            ks = wavenumbers(ωs[i], micro;
                symmetry = symmetry,
                num_wavenumbers = num_wavenumbers,
                bisection_method = bisection_method,
                tol = tol,
                k_effs = Complex{T}[keffs[i-1]],
            kws...)

            if isempty(ks)
                error("No effective wavenumbers were found for the frequency ω = $(ωs[i]), the $(i)th frequency given, so the wavenumber can not be tracked any further.")
            end

            j = findmin(norm.(ks .- keffs[i-1]))[2]

            if norm(ks[j] - keffs[i-1]) > kscale
                @warn "The wavenumber found for the frequency ω = $(ωs[i]), the $(i)th frequency given, is a distance $(norm(ks[j] - keffs[i-1])) from the previous wavenumber, which is larger than the background wavenumber $(kscale). The tracking may have jumped to a different branch."
            end

            keffs[i] = ks[j]
            errors[i] = disp(keffs[i]) |> abs

        elseif interpolation_order == 1

            k1, k2 = keffs[i-2], keffs[i-1]

            # Linear interpolation: k3_predicted = k1 + (k2 - k1) * (ω3 - ω1) / (ω2 - ω1)
            keffs[i] = k1 + (k2 - k1) * (ωs[i] - ωs[i-2]) / (ωs[i-1] - ωs[i-2])

            refine!(i, disp, kscale)

        elseif interpolation_order == 2
            # Quadratic interpolation using the three previous nodes: keffs[i-3], keffs[i-2], keffs[i-1]
            # Fit a quadratic polynomial through (ωs[i-3], keffs[i-3]), (ωs[i-2], keffs[i-2]), (ωs[i-1], keffs[i-1])
            # and evaluate at ωs[i]

            ω1, ω2, ω3 = ωs[i-3], ωs[i-2], ωs[i-1]
            k1, k2, k3 = keffs[i-3], keffs[i-2], keffs[i-1]
            ω_target = ωs[i]

            # Lagrange interpolation formula for quadratic interpolation
            L1 = ((ω_target - ω2) * (ω_target - ω3)) / ((ω1 - ω2) * (ω1 - ω3))
            L2 = ((ω_target - ω1) * (ω_target - ω3)) / ((ω2 - ω1) * (ω2 - ω3))
            L3 = ((ω_target - ω1) * (ω_target - ω2)) / ((ω3 - ω1) * (ω3 - ω2))

            keffs[i] = L1 * k1 + L2 * k2 + L3 * k3

            refine!(i, disp, kscale)

        else
            throw(ArgumentError("interpolation_order = $(interpolation_order) is not supported, it needs to be 0, 1, or 2."))
        end

        return keffs[i], errors[i]
    end

    interpolation_order = 0
    for i in 2:length(ωs)

        # an interpolation of order n needs the n + 1 previous wavenumbers
        interpolation_order = min(interpolation_order, i - 2)

        predict_next!(interpolation_order, i; basis_order = basis_orders[i], kws...)

        k = ωs[i] / real(micro.medium.c)

        interpolation_order = if abs(keffs[i] - keffs[i-1]) / k  < 1.0
            min(interpolation_order + 1, 2)
        elseif abs(keffs[i] - keffs[i-1]) / k  < 3.0
            1
        else 0
        end

        interpolation_orders[i] = interpolation_order
    end

    @debug "The largest residual abs(det(MM(k_eff))) over all the frequencies was $(maximum(errors)), and the interpolation orders used were $(interpolation_orders)."

    keffs = map(keffs) do keff
        imag(keff) < -tol ? - keff : keff
    end

    # kφs = map(ωs) do ω
    #     wavenumber_low_volumefraction(ω, micro; basis_order = basis_orders[i],  kws...)
    # end

    # plot(ωs, real.(keffs), lab = "real")
    # plot!(ωs, real.(kφs), lab = "real (low volume fraction)", linestyle=:dash)
    # plot(ωs, imag.(keffs), lab = "imag")
    # plot!(ωs, imag.(kφs), lab = "imag (low volume fraction)", linestyle=:dash)
    # plot!(ωs, interpolation_orders, lab = "order")
    # plot!(ωs, errors / tol, lab = "error")

    return keffs
end
