import StaticArrays: SVector
import EffectiveWaves: AbstractAzimuthalSymmetry

using HCubature, Interpolations, Statistics

function discrete_system_residue(discrete_coefs, ω::T, source::AbstractSource{Acoustic{T,Dim}}, material::Material{Sphere{T,Dim}}, ::AbstractAzimuthalSymmetry{Dim};
        basis_order::Int = 1,
        mesh_points::Int = 11,
        rtol::T = 1e-2,
        maxevals::Int = Int(1e5),
        pair_corr = hole_correction_pair_correlation
    ) where {T,Dim}

    if length(material.microstructure.species) > 1
        @warn "discrete_system has only been implemented for 1 species for now. Will use only first specie."
    end

    s1 = material.microstructure.species[1]

    R = outer_radius(material.shape)

    gs = regular_spherical_coefficients(source)(basis_order,origin(material.shape),ω);
    v = regular_basis_function(source.medium,  ω)

    Uout = outgoing_translation_matrix(ω, material;
        basis_order = basis_order, tol = rtol
    )

    t_matrices = get_t_matrices(source.medium, material.microstructure.species, ω, basis_order)
    t_diags = diag.(t_matrices)

    rθφ2xyz = radial_to_cartesian_coordinates

    len = basisorder_to_basislength(Acoustic{T,Dim}, basis_order)

    function incident_coefficients(rθ1::AbstractVector{T})
        lm2n = lm_to_spherical_harmonic_index
        vs = v(2basis_order, rθφ2xyz([rθ1; 0.0]))

        coefs = [
            t_diags[1][lm2n(l,m)] *
            sum(
                gaunt_coefficient(dl,0,l,m,l1,-m) * vs[lm2n(l1,-m)] * gs[lm2n(dl,0)]
            for l1 in abs(m):(2basis_order), dl in 0:basis_order)
        for l = 0:basis_order for m = -l:l]

        return coefs
    end

    ls, ms = spherical_harmonics_indices(basis_order)

    function kernel_function(rθ1::AbstractVector{T})
        x1 = rθφ2xyz(SVector(rθ1[1],rθ1[2],zero(T)))

        fun = function (rθφ::SVector{3,T})
            x2 = rθφ2xyz(rθφ)
            if pair_corr(x1,s1,x2,s1) ≈ zero(T)
                return zeros(Complex{T}, len)
            end
            # U = outgoing_translation_matrix(medium, basis_order, basis_order, ω, x1 - x2)
            U = Uout(x1 - x2)

            U = U .* (number_density(s1) * pair_corr(x1,s1,x2,s1) * sin(rθφ[2]) * rθφ[1]^2)
            return [
                    t_diags[1][n] * sum(U[nd,n] * discrete_coefs(x2)[nd] for nd in 1:len)
            for n in 1:len]
        end

        return fun
    end

    rs = LinRange(0.0, R - outer_radius(s1), mesh_points);
    θs = LinRange(0.0, π, mesh_points);

    test_ker = kernel_function(SVector(mean(rs),θs[1]))
    (v2,err) = hcubature(test_ker, SVector(0.0,0.0,-π), SVector(R-outer_radius(s1),π,π)
        ; rtol=rtol, maxevals = maxevals
    );

    println("The estimated max coefficient of the integrated kernel is:")
    println("I = ", maximum(abs.(v2)) )
    println("with an estimated error of: ")
    println("E = ", err)


    relative_residues = [
        begin
            ker = kernel_function([r,θ])
            K = hcubature(ker, SVector(0.0,0.0,-π), SVector(R-outer_radius(s1),π,π)
                ; rtol=rtol, maxevals = maxevals
            )[1];
            inc = incident_coefficients([r,θ]);
            norm(- discrete_coefs(rθφ2xyz([r,θ,0.0])) + inc + K) / norm(inc)
        end
    for r in rs, θ in θs]

    return (mean(relative_residues),maximum(relative_residues))
end
