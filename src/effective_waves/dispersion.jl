# eigensystem(ω::AbstractFloat, medium::PhysicalMedium, species::Species, sym::AbstractSymmetry; kws...) = eigensystem(ω, Microstructure(medium,species), sym; kws...)

dispersion_equation(ω::AbstractFloat, source::AbstractSource, material::Material; kws...) = dispersion_equation(ω, material.microstructure, Symmetry(source,material); kws...)

dispersion_equation(ω, medium::PhysicalMedium, sps::Species, symmetry::AbstractSymmetry; kws...) = dispersion_equation(ω, Microstructure(medium,sps), symmetry; kws...)

dispersion_complex(ω, medium::PhysicalMedium, sps::Species, symmetry::AbstractSymmetry; kws...) = dispersion_complex(ω, Microstructure(medium,sps), symmetry; kws...)

function dispersion_equation(ω::T, micro::Microstructure{Dim}, symmetry::AbstractSymmetry = PlanarSymmetry{Dim}();
        basis_order = 3 * Int(round(maximum(outer_radius.(micro.species)) * ω / abs(micro.medium.c) )) + 1,
        tol::T = 1e-4, low_tol::T = max(1e-4, tol), kws...
    ) where {T<:Number, Dim}

    # low_tol: a tolerance used for a first pass with time_limit

    MM = eigensystem(ω, micro, symmetry; basis_order=basis_order, kws... )

    # the constraint uses keff_vec[2] < -low_tol to better specify solutions where imag(k_effs)~0 and imag(k_effs)<0
    constraint(k_eff::Complex{T}) = (imag(k_eff) < -low_tol) ? (-one(T) + exp(-T(100.0) * imag(k_eff))) : zero(T)

    function detMM(k_eff::Complex{T})
        constraint(k_eff) + abs2(det(MM(k_eff)))
    end

    return detMM
end

function dispersion_complex(ω::T, micro::Microstructure{Dim}, symmetry::AbstractSymmetry = PlanarSymmetry{Dim}();
        basis_order = 3 * Int(round(maximum(outer_radius.(micro.species)) * ω / abs(micro.medium.c) )) + 1,
        tol::T = 1e-4, kws...
    ) where {T<:Number, Dim}

    MM = eigensystem(ω, micro, symmetry; basis_order=basis_order, kws... )
    detMM(k_eff::Complex{T})::Complex{T} = det(MM(k_eff))

    return detMM
end
