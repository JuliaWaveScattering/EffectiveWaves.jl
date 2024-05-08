
function material_scattering_coefficients(scat_field::ScatteringCoefficientsField;
        numdensity = (x1, s1) -> number_density(s1),
        rtol::AbstractFloat = 1e-2,
        maxevals::Int = Int(5e4)
    )

    lm2n = lm_to_spherical_harmonic_index
    rθφ2xyz = radial_to_cartesian_coordinates

    if scat_field.medium != scat_field.material.microstructure.medium @error mismatched_medium end

    v = regular_basis_function(scat_field.medium, scat_field.ω)

    particle_radius = maximum(outer_radius.(scat_field.material.microstructure.species))
    R = outer_radius(scat_field.material.shape)

    function kernel(rθφ)
        x = rθφ2xyz(rθφ)

        vs = conj.(v(scat_field.basis_order + scat_field.basis_field_order, x))
        fs = scat_field.coefficient_field(x)

        (rθφ[1]^2 * sin(rθφ[2]) * numdensity(x,scat_field.material.microstructure.species[1])) .* [
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
        numdensity = (x1, s1) -> number_density(s1),
        rtol::AbstractFloat = 1e-4,
        maxevals::Int = Int(1e5),
    ) where T

    if scat_field.medium != scat_field.material.microstructure.medium @error mismatched_medium end

    k = scat_field.ω / scat_field.medium.c
    R = outer_radius(scat_field.material.shape)

    lm2n = lm_to_spherical_harmonic_index
    ls = 0:scat_field.basis_order

    rθφ2xyz = radial_to_cartesian_coordinates

    particle_radius = maximum(outer_radius.(scat_field.material.microstructure.species))

    function kernel(r)
        x = rθφ2xyz([r,0.0,0.0])

        js = sbesselj.(ls, k*r)

        fs = scat_field.coefficient_field(x)

        (4.0 * π * r^2 * numdensity([r,0,0],scat_field.material.microstructure.species[1])) .* sum(
            (-1.0) .^ ls .* sqrt.(2.0 .* ls .+ 1.0) .* fs[lm2n.(ls,0)] .* js
        )
    end

    (v,err) = hquadrature(kernel, 0.0, R - particle_radius; rtol=rtol, maxevals = maxevals)

    println("Integration error of $err")

    return [v]

end

function material_effective_tmatrix(ω::T, material::Material{Sphere{T,2}};
    basis_field_order::Int = 2,
    rtol::AbstractFloat = 1e-4,
    basis_order::Int = 2basis_field_order
    ) where T

    # Extracting important parameters
    medium = material.microstructure.medium
    k = ω / medium.c
    R = outer_radius(material.shape)
    species = material.microstructure.species
    particle_radius = maximum(outer_radius.(species))
    Rtildas = R .- outer_radius.(species)

    # Computing effective wavenumber
    k_eff = wavenumbers(ω, medium, species; basis_order = basis_order, num_wavenumbers = 5)[1]

    # Solving artificial eigensystem
    F = eigenvectors(ω, k_eff, material.microstructure, PlanarSymmetry{2}(); basis_order = basis_order)

    nbo, n_λ, _ = size(F)
    basis_order = Int((nbo-1)/2)
    L = basis_order+basis_field_order

    n_densities = [number_density(sp) for sp in species]

    J = [besselj(n,k*Rtildas[λ])*n_densities[λ] 
        for n = -L-1:L, λ in 1:n_λ]

    Jstar = [besselj(n,k_eff*Rtildas[λ])*n_densities[λ] 
        for n = -L-1:L, λ in 1:n_λ]

    H = [hankelh1(n,k*Rtildas[λ])*n_densities[λ] 
        for n = -L-1:L, λ in 1:n_λ]

    pre_num = k*J[1:1+2*L,:].*Jstar[2:2*(1+L),:] - k_eff*J[2:2*(1+L),:].*Jstar[1:1+2*L,:]
    pre_denom = k*H[1:1+2*L,:].*Jstar[2:2*(1+L),:] - k_eff*H[2:2*(1+L),:].*Jstar[1:1+2*L,:]

    Matrix_Num = complex(zeros(2*basis_order+1,2*basis_field_order+1))
    Matrix_Denom = complex(zeros(2*basis_order+1,2*basis_field_order+1))
    for n = 1:1+2*basis_field_order        
        Matrix_Num[:,n] = vec(sum(pre_num[n+2*basis_order:-1:n,:].*F,dims=2))
        Matrix_Denom[:,n] = vec(sum(pre_denom[n+2*basis_order:-1:n,:].*F,dims=2))
    end

    Tmats = (- vec(sum(Matrix_Num,dims=1)./sum(Matrix_Denom,dims=1)))

    return Tmats

end
