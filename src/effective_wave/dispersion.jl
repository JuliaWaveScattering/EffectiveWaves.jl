dispersion_equation(ω::AbstractFloat, source::AbstractSource, material::Material; kws...) = dispersion_equation(ω, source.medium, material.species, Symmetry(source,material); kws...)


function dispersion_equation(ω::T, medium::PhysicalMedium{Dim}, species::Species{Dim}, symmetry::AbstractSymmetry = PlanarSymmetry{Dim}();
        basis_order = 3 * Int(round(maximum(outer_radius.(species)) * ω / abs(medium.c) )) + 1,
        tol::T = 1e-4, low_tol::T = max(1e-4, tol), kws...
    ) where {T<:Number, Dim}

    # low_tol: a tolerance used for a first pass with time_limit

    MM = eigensystem(ω, medium, species, symmetry; basis_order=basis_order, kws... )

    # the constraint uses keff_vec[2] < -low_tol to better specify solutions where imag(k_effs)~0 and imag(k_effs)<0
    constraint(k_eff::Complex{T}) = (imag(k_eff) < -low_tol) ? (-one(T) + exp(-T(100.0) * imag(k_eff))) : zero(T)

    function detMM(k_eff::Complex{T})
        constraint(k_eff) + abs(det(MM(k_eff)))
    end

    return detMM
end

function dispersion_complex(ω::T, medium::PhysicalMedium{Dim}, species::Species{Dim}, symmetry::AbstractSymmetry = PlanarSymmetry{Dim}();
        basis_order = 3 * Int(round(maximum(outer_radius.(species)) * ω / abs(medium.c) )) + 1,
        tol::T = 1e-4, kws...
    ) where {T<:Number, Dim}

    MM = eigensystem(ω, medium, species, symmetry; basis_order=basis_order, kws... )
    detMM(k_eff::Complex{T})::Complex{T} = det(MM(k_eff))

    return detMM
end
