# NOTE: PlanarAzimuthalSymmetry() does not include all possible wavenumbers
function wavenumbers(ωs::AbstractVector{T}, micro::Microstructure{Dim};
        tol::T = 1e-5, branch_number::Int = 1, kws...) where {T,Dim}

    a = mean(micro.species .|> outer_radius)
    kas = a .* ωs ./ micro.medium.c

    if abs(kas[1]) > 0.1
        @warn "The wavenumber is large compared to the particle size with ka = $(kas[1]). Consider using a smaller frequency for this method which sweeps over frequencies."
    end

    k0s = wavenumbers_path(ωs[1], micro; num_wavenumbers = 1, tol = tol, kws...)

    if branch_number > length(k0s)
        @error "Branch number $(branch_number) is larger than the number of wavenumbers found $(length(k0s)). Returning the last wavenumber found."
        branch_number = length(k0s)
    end
    keffs = Vector{Complex{T}}(undef, length(ωs))
    keffs[1] = k0s[branch_number]
    

end