dispersion_equation(ω::AbstractFloat, source::AbstractSource, material::Material; kws...) = dispersion_equation(ω, source.medium, material.species, setupsymmetry(source,material); kws...)


function dispersion_equation(ω::T, medium::PhysicalMedium{T,Dim}, species::Species{T,Dim}, symmetry::AbstractSetupSymmetry = PlanarSymmetry(); tol::T = 1e-4, kws...) where {T<:Number, Dim}

    low_tol = max(1e-4, tol) # a tolerance used for a first pass with time_limit

    MM = eigensystem(ω, medium, species, symmetry; kws... )

    # the constraint uses keff_vec[2] < -low_tol to better specify solutions where imag(k_effs)~0 and imag(k_effs)<0
    constraint(keff_vec::Vector{T}) = (keff_vec[2] < -low_tol) ? (-one(T) + exp(-T(100.0)*keff_vec[2])) : zero(T)

    function detMM(keff_vec::Vector{T})
        constraint(keff_vec) + sqrt(abs(det(MM(keff_vec[1]+im*keff_vec[2]))))
    end

    return detMM
end
