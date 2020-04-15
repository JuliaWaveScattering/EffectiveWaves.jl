# Here we calculate the effective wavenumber and effective wave amplitudes, without any restriction to the volume fraction of particles and incident wave frequency.

# include depricated function to find a single effective wavenumber, when in fact there are many. The code is still used in tests and gives many correct results
# include("wavenumber_single.jl")

" Returns all the transmitted effective wavenumbers"
wavenumbers(ω::T, medium::PhysicalMedium{T}, specie::Specie{T}; kws...) where T<:Number = wavenumbers(ω, medium, [specie]; kws...)

function wavenumbers(ω::T, medium::PhysicalMedium{T}, species::Species{T};
        num_wavenumbers::Int = 2, tol::T = 1e-5,
        kws...) where T<:Number

    # for very low attenuation, need to search close to assymptotic root with a path method.

    k_effs = wavenumbers_path(ω, medium, species;
    num_wavenumbers = 2, tol = tol, kws...)

    if num_wavenumbers > 2
        k_effs2 = wavenumbers_bisection(ω, medium, species;
            num_wavenumbers=num_wavenumbers, tol = tol, kws...)
        k_effs = [k_effs; k_effs2]
        k_effs = reduce_kvecs(k_effs, tol)
        k_effs = sort(k_effs, by = imag)
    end

    return k_effs
end
