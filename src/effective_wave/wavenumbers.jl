# Here we calculate the effective wavenumber and effective wave amplitudes, without any restriction to the volume fraction of particles and incident wave frequency.

" Returns all the transmitted effective wavenumbers"
wavenumbers(ω::T, medium::PhysicalMedium{T}, specie::Specie{T}; kws...) where T<:Number = wavenumbers(ω, medium, [specie]; kws...)

function wavenumbers(ω::T, medium::PhysicalMedium{T}, species::Species{T};
        num_wavenumbers::Int = 2, tol::T = 1e-5,
        max_Imk::T = zero(T), max_Rek::T = zero(T),
        kws...) where T<:Number

    # For very low attenuation, need to search close to assymptotic root with a path method.
    k_effs = wavenumbers_path(ω, medium, species;
    num_wavenumbers = 2, max_Imk=max_Imk, tol = tol, kws...)

    # NOTE: these search methods would significantly improve if we used the asymptotic result for multiple wavenumbers and monopole scatterers. This would give a reasonable length scale and on where to search.
    if num_wavenumbers > 2
        box_k = box_keff(ω, medium, species; tol = tol)
        max_imag = max(3.0 * maximum(imag.(k_effs)), max_Imk)
        max_imag = max(max_imag, box_k[2][2])
        max_real = max(2.0 * maximum(real.(k_effs)), max_Rek)
        max_real = max(max_real, box_k[1][2])
        box_k = [[-max_real,max_real], [0.0,max_imag]]

        k_effs2 = wavenumbers_bisection(ω, medium, species;
            # num_wavenumbers=num_wavenumbers,
            tol = tol, box_k = box_k,
            kws...)
        k_effs = [k_effs; k_effs2]
        k_effs = reduce_kvecs(k_effs, tol)
        k_effs = sort(k_effs, by = imag)
    end

    return k_effs
end
