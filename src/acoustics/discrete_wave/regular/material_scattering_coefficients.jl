
function material_scattering_coefficients(scat_field::ScatteringCoefficientsField;
        rtol::AbstractFloat = 1e-2,
        maxevals::Int = Int(5e4)
    )

    lm2n = lm_to_spherical_harmonic_index
    rθφ2xyz = radial_to_cartesian_coordinates

    v = regular_basis_function(scat_field.medium, scat_field.ω)

    numdensity = number_density(scat_field.material.species)

    particle_radius = maximum(outer_radius.(scat_field.material.species))

    function kernel(rθφ)
        x = rθφ2xyz(rθφ)

        vs = conj.(v(scat_field.basis_order + scat_field.basis_field_order, x))
        fs = scat_field.coefficient_field(x)

        (rθφ[1]^2 * sin(rθφ[2]) * numdensity) .* [
            sum(
                gaunt_coefficient(l,m,dl,dm,l1,m-dm) * vs[lm2n(l1,m-dm)] * fs[lm2n(dl,dm)]
            for dl = 0:scat_field.basis_order for dm = -dl:dl for l1 in max(abs(m-dm),abs(dl-l)):(dl+l))
        for l = 0:scat_field.basis_field_order for m = -l:l]
    end

    (v,err) = hcubature(kernel, SVector(0.0,0.0,-π), SVector(R - particle_radius,π,π)
        ; rtol=rtol, maxevals = maxevals
    )

    println("Integration error of $err")

    return v
end

function material_scattering_coefficients(scat_field::ScatteringCoefficientsField{Sphere{T,3},Acoustic{T,3},RadialSymmetry{3}};
        rtol::AbstractFloat = 1e-4,
        maxevals::Int = Int(1e5)
    ) where T

    k = scat_field.ω / scat_field.medium.c
    R = outer_radius(scat_field.material.shape)

    lm2n = lm_to_spherical_harmonic_index
    ls = 0:scat_field.basis_order

    rθφ2xyz = radial_to_cartesian_coordinates

    numdensity = number_density(scat_field.material.species)
    particle_radius = maximum(outer_radius.(scat_field.material.species))

    function kernel(r)
        x = rθφ2xyz([r,0.0,0.0])

        js = sbesselj.(ls, k*r)

        fs = scat_field.coefficient_field(x)

        (4.0 * π * r^2 * numdensity) .* sum(
            (-1.0) .^ ls .* sqrt.(2.0 .* ls .+ 1.0) .* fs[lm2n.(ls,0)] .* js
        )
    end

    (v,err) = hquadrature(kernel, 0.0, R - particle_radius; rtol=rtol, maxevals = maxevals)

    println("Integration error of $err")

    return [v]

end
