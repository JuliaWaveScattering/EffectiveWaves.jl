"""
    wavenumbers(ωs::AbstractVector{T}, micro::Microstructure{Dim}; symmetry, tol, branch_number, optimoptions, kws...)

Compute effective wavenumbers over a range of frequencies using adaptive interpolation and refinement.

## Arguments
- `tol::T`: Tolerance for numerical optimization (default: `1e-6`)
- `branch_number::Int`: Which branch of wavenumbers to return (default: `1`)

## Returns
- `keffs::Vector{Complex{T}}`: Complex effective wavenumbers with positive imaginary part

## Method
- The first effective wavenumber is computed at ω₁ using exact root-finding
- The second wavenumber is computed at ω₂ using exact root-finding
- `keffs[3]` is predicted using linear interpolation of `keffs[1]` and `keffs[2]`
- For i > 3, `keffs[i]` is predicted using quadratic (Lagrange) interpolation of the three previous nodes
- Each prediction is refined by optimization to satisfy the dispersion relation

## Notes
- The method assumes small step sizes in wavenumber space (Δka < 0.001)
- Starting frequencies should be small compared to particle size (ka < 0.1)
- Physical solutions must have Im(k_eff) > 0; negative imaginary parts are negated
- NOTE: PlanarAzimuthalSymmetry() does not include all possible wavenumbers
"""
function wavenumbers(ωs::AbstractVector{T}, micro::Microstructure{Dim}; 
        symmetry::AbstractSymmetry{Dim} = PlanarAzimuthalSymmetry{Dim}(), 
        tol::T = 1e-6, branch_number::Int = 1,
        basis_orders = 2 * ones(Int, length(ωs)),
        optimoptions = Optim.Options(
            iterations = Int(round(-log(tol))) * 20,
            g_tol = tol^T(2), x_abstol=tol^T(2)),
        kws...) where {T,Dim}

    a = mean(micro.species .|> outer_radius)
    kas = a .* ωs ./ real(micro.medium.c)
    dka = abs.((kas - circshift(kas,1))[2:end]) |> minimum

    if dka > 0.001
        @warn "The method is design to take very small step increases in the wavenumber ka. The minimum step increase in ka is $(dka), but should be smaller than 0.001."
    end

    if abs(kas[1]) > 0.1
        @warn "The method is design to start with small wavelengths, compared to the particle size. Starting with ka = $(kas[1]) might be too large. "
    end

    function predict_next!(interpolation_order::Int, i::Int, ωs::AbstractVector{T}, keffs::Vector{Complex{T}}, errors::Vector{T}, micro::Microstructure{Dim}, symmetry::AbstractSymmetry{Dim}, tol::T; kws...) where {T,Dim}
        
        disp = dispersion_complex(ωs[i], micro, symmetry; kws...)
        f_vec(x_vec) = abs(disp(x_vec[1] + x_vec[2]*im))

        if interpolation_order == 0
            ks = wavenumbers(ωs[i], micro; 
                symmetry = symmetry,
                num_wavenumbers = 1, 
                tol = tol, 
                k_effs = Complex{T}[keffs[i-1]], 
            kws...)
            j = findmin(norm.(ks .- keffs[i-1]))[2]

            keffs[i] = ks[j]
            errors[i] = disp(keffs[i]) |> abs

        elseif interpolation_order == 1

            k1, k2 = keffs[i-2], keffs[i-1]

            # Linear interpolation: k3_predicted = k1 + (k2 - k1) * (ω3 - ω1) / (ω2 - ω1)
            keffs[i] = k1 + (k2 - k1) * (kas[i] - kas[i-2]) / (kas[i-1] - kas[i-2])

            res = optimize(f_vec, [keffs[i] |> real, keffs[i] |> imag], NelderMead(), optimoptions);

            keffs[i] = res.minimizer[1] + res.minimizer[2]*im
            errors[i] = disp(keffs[i]) |> abs

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

            res = optimize(f_vec, [keffs[i] |> real, keffs[i] |> imag], NelderMead(), optimoptions);
            keffs[i] = res.minimizer[1] + res.minimizer[2]*im
            errors[i] = disp(keffs[i]) |> abs
        end

        return keffs[i], errors[i]
    end

    # ceff = effective_medium(micro).c
    disp = dispersion_complex(ωs[1], micro, symmetry; basis_order = basis_orders[1], kws...)
    f_vec(x_vec) = abs(disp(x_vec[1] + x_vec[2]*im))

    k0s = wavenumbers(ωs[1], micro; symmetry = symmetry, basis_order = basis_orders[1], num_wavenumbers = 1, tol = tol, kws...)
    # k0s = wavenumbers_bisection_robust(ωs[1], micro; symmetry = symmetry, num_wavenumbers = 10, tol = tol, kws...)
    
    if branch_number > length(k0s)
        @error "Branch number $(branch_number) is larger than the number of wavenumbers found $(length(k0s)). Returning the last wavenumber found."
        branch_number = length(k0s)
    end

    errors = Vector{T}(undef, length(ωs))
    
    keffs = Vector{Complex{T}}(undef, length(ωs))
    keffs[1] = k0s[branch_number]
    errors[1] = disp(keffs[1]) |> abs

    interpolation_orders = zeros(Int, length(ωs))
    interpolation_orders[1] = 0

    interpolation_order = 0
    for i in 2:length(ωs)

        predict_next!(interpolation_order, i, ωs, keffs, errors, micro, symmetry, tol;                 basis_order = basis_orders[i], kws...)
        keffs[i]
        k = ωs[i] / real(micro.medium.c)
        abs(keffs[i] - keffs[i-1]) / k 

        interpolation_order = if abs(keffs[i] - keffs[i-1]) / k  < 1.0 
            min(interpolation_order + 1, 2)
        elseif abs(keffs[i] - keffs[i-1]) / k  < 3.0
            1
        else 0
        end

        interpolation_orders[i] = interpolation_order
    end

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