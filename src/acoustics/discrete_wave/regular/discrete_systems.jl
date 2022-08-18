# using EffectiveWaves, LinearAlgebra
# import StaticArrays: SVector

import MultipleScattering: outgoing_translation_matrix

function discrete_system(ω::T, source::AbstractSource{Acoustic{T,Dim}}, material::Material{Dim,Sphere{T,Dim}}, ::WithoutSymmetry{Dim}; kws...) where {T,Dim}

    return discrete_system(ω, source, material, AzimuthalSymmetry{Dim}(); kws...)
end

"""
    discrete_system(k,r1_vec::Vector,r2,θ2; basis_order::Int = 3)

documentation
"""
function discrete_system(ω::T, source::AbstractSource{Acoustic{T,Dim}}, material::Material{Dim,Sphere{T,Dim}}, ::AbstractAzimuthalSymmetry{Dim};
        basis_order::Int = 1,
        basis_field_order::Int = 2,
        legendre_order::Int = basis_field_order + 1,
        mesh_points::Int = Int(round(sqrt(1.1 * (basis_field_order) * legendre_order ))) + 2,
        rtol::T = 1e-2,
        maxevals::Int = Int(2e4),
        numdensity = (x1, s1) -> number_density(s1),
        pair_corr = hole_correction_pair_correlation
    ) where {T,Dim}

    if length(material.species) > 1
        @warn "discrete_system has only been implemented for 1 species for now. Will use only first specie."
    end

    s1 = material.species[1]
    scale_number_density = one(T) - one(T) / material.numberofparticles
    bar_numdensity(x1,s1) = scale_number_density * numdensity(x1,s1)

    R = outer_radius(material.shape)

    gs = regular_spherical_coefficients(source)(basis_field_order,origin(material.shape),ω);

    v = regular_basis_function(source.medium,  ω)

    Uout = outgoing_translation_matrix(ω, source.medium, material;
        basis_order = basis_order, tol = rtol
    )

    t_matrices = get_t_matrices(source.medium, material.species, ω, basis_order)
    t_diags = diag.(t_matrices)

    rθφ2xyz = radial_to_cartesian_coordinates

    r1s = LinRange(0,R-outer_radius(s1), mesh_points)
    θ1s = LinRange(0,π, mesh_points)

    ls, ms = spherical_harmonics_indices(basis_order)
    azi_inds(m::Int) = lm_to_spherical_harmonic_index.(abs(m):basis_field_order,-m)

    len = basisorder_to_basislength(Acoustic{T,Dim}, basis_order)
    len_p = sum(legendre_order for nd in 1:len for i in azi_inds(ms[nd]))

    function incident_coefficients(r1s::AbstractVector{T},θ1s::AbstractVector{T})
        lm2n = lm_to_spherical_harmonic_index

        coefs = [
            begin
                vs = v(basis_order + basis_field_order, rθφ2xyz(SVector(r1,θ1,zero(T))))
                data = [
                    t_diags[1][lm2n(l,m)] *
                    sum(
                        gaunt_coefficient(dl,0,l,m,l1,-m) * vs[lm2n(l1,-m)] * gs[lm2n(dl,0)]
                    for dl in 0:basis_field_order for l1 in max(abs(m),abs(dl-l)):(dl+l))
                for l = 0:basis_order for m = -l:l]

                data
            end
        for r1 in r1s, θ1 in θ1s][:];

        return vcat(coefs...)
    end

    function field_basis(rθφ::AbstractVector{T})
        Ys = spherical_harmonics(basis_field_order, rθφ[2], rθφ[3]);
        P = Legendre{T}()

        # [P_0(cos(θ)), …, P_(legendre_order-1)(cos(θ))]
        Prs = P[2 * rθφ[1] / (R - outer_radius(s1)) - one(T), 1:legendre_order]

        return [Pr * Yθ for Pr in Prs, Yθ in Ys]
    end
    # function field_basis(rθ::AbstractVector{T})
    #     P = Legendre{T}()
    #
    #     # [P_0(cos(θ)), …, P_(legendre_order-1)(cos(θ))]
    #     Pθs = P[cos(rθ[2]), 1:legendre_order]
    #     Prs = P[2 * rθ[1] / (R - outer_radius(s1)) - one(T), 1:legendre_order]
    #
    #     # [Pr * Pθ for Pr in Prs, Pθ in Pθs][:]
    #     return (Prs * transpose(Pθs))[:]
    # end

    function kernel_function(rθ1::SVector{2,T})
        x1 = rθφ2xyz(SVector(rθ1[1],rθ1[2],zero(T)))
        Kzero = zeros(Complex{T},len,len_p)

        fun = function (rθφ::SVector{3,T})
            x2 = rθφ2xyz(rθφ)
            if pair_corr(x1,s1,x2,s1) ≈ zero(T)
                return Kzero
            end
            basis2 = field_basis(rθφ)

            U = Uout(x1 - x2)
            U = U .* (bar_numdensity(x2,s1) * pair_corr(x1,s1,x2,s1) * sin(rθφ[2]) * rθφ[1]^2)

            return reshape(
                [
                    t_diags[1][n] * U[nd,n] * b2
                for nd in 1:len for b2 in basis2[:,azi_inds(ms[nd])][:] for n in 1:len],
            (len, :))
        end

        return fun
    end

    function δφj(rθφ1::AbstractVector{T})
        basis1 = field_basis(rθφ1)
        return reshape(
            [
                (nd == n) ? b1 : zero(Complex{T})
            for nd in 1:len for b1 in basis1[:,azi_inds(ms[nd])][:] for n in 1:len],
        (len, len_p))
    end

    δφj(rθ1::SVector{2,T}) = δφj(SVector(rθ1[1],rθ1[2],zero(T)))

    test_ker = kernel_function(SVector(mean(r1s),θ1s[1]))
    (I,E) = hcubature(test_ker, SVector(0.0,0.0,-π), SVector(R-outer_radius(s1),π,π);
        rtol=rtol, maxevals=maxevals
    );

    println("The estimated max coefficient of the integrated kernel is:")
    println("I = ", maximum(abs.(I)) )
    println("with an estimated error of: ")
    println("E = ", E)

    Ks = [
        begin
            rθ1 = SVector(r1,θ1)
            ker = kernel_function(rθ1)
            ker_integrated = hcubature(ker, SVector(0.0,0.0,-π), SVector(R-outer_radius(s1),π,π);
                rtol=rtol, maxevals=maxevals
            )[1]

            δφj(rθ1) - ker_integrated
        end
    for r1 in r1s, θ1 in θ1s][:];

    bigK = vcat(Ks...);

    bs = incident_coefficients(r1s,θ1s);

    Fs = bigK \ bs;

    ## Alternative:
    # as = inv(transpose(conj.(bigK)) * bigK) * transpose(conj.(bigK)) * bs;

    if norm(bigK * Fs - bs) / norm(bs) > rtol
        @warn "Numerical solution has a relative residual error of $(norm(bigK * Fs - bs) / norm(bs)), where the requested relative tolernance was: $rtol. This residual error can be decreased by increasing the basis_field_order (current value: $basis_field_order) for the field."
    end

    function scattered_field(xs::AbstractVector{T})
        rθφ = cartesian_to_radial_coordinates(xs)
        return δφj(rθφ) * Fs
        # The factor exp(-im * m * φ) is due to azimuthal symmetry
        # azi_factor = exp.((-im*rθφ[3]) .* ms)
        # return azi_factor .* (as * field_basis(rθφ[1:2]))
    end

    return ScatteringCoefficientsField(ω, source.medium, material, scattered_field;
        symmetry = AzimuthalSymmetry{Dim}(),
        basis_order = basis_order,
        basis_field_order = basis_field_order
    )
end

function discrete_system(ω::T, source::AbstractSource{Acoustic{T,3}}, material::Material{3,Sphere{T,3}}, ::RadialSymmetry{3};
        basis_order::Int = 1,
        basis_field_order::Int = Int(round(T(2) * real(ω / source.medium.c) * outer_radius(material.shape))) + 1,
        legendre_order::Int = 2basis_field_order + 1,
        # mesh_points::Int = 3 * (legendre_order + 1),
        mesh_points::Int = 3 * (legendre_order + 1),
        rtol::T = 1e-2, ptol::T = 1e-2,
        maxevals::Int = Int(2e4),
        numdensity = (x1, s1) -> number_density(s1),
        pair_corr = hole_correction_pair_correlation
    ) where T

    if length(material.species) > 1
        @warn "discrete_system has only been implemented for 1 species for now. Will use only first specie."
    end

    s1 = material.species[1]
    a12 = 2.0 * s1.exclusion_distance * outer_radius(s1);

    scale_number_density = one(T) - one(T) / material.numberofparticles
    bar_numdensity(x1,s1) = scale_number_density * numdensity(x1,s1)

    R = outer_radius(material.shape)

    g0 = regular_spherical_coefficients(source)(1,origin(material.shape),ω)[1];

    v = regular_basis_function(source.medium,  ω)

    Uout = outgoing_translation_matrix(ω, source.medium, material;
        basis_order = basis_order, tol = rtol
    )

    lm_to_n = lm_to_spherical_harmonic_index
    ls_to_ns = lm_to_n.(0:basis_order,0)

    t_matrices = get_t_matrices(source.medium, material.species, ω, basis_order)
    t_diags = diag.(t_matrices)

    rθφ2xyz = radial_to_cartesian_coordinates

    r1s = LinRange(0,R-outer_radius(s1), mesh_points)

    len = basis_order + 1
    len_p = legendre_order * len

    function incident_coefficients(r1s::AbstractVector{T})

        coefs = [
            begin
                vs = v(basis_order, rθφ2xyz([r1,zero(T),zero(T)]))
                t_diags[1][ls_to_ns] .* vs[ls_to_ns] .* (-T(1)) .^ (0:basis_order)
            end
        for r1 in r1s][:];

        return vcat((sqrt(4pi) * g0) .* coefs...)
    end

    function field_basis(rθφ::AbstractVector{T})
        Ys = conj.(spherical_harmonics(basis_order, rθφ[2], rθφ[3]));
        P = Legendre{T}()

        Prs = P[2 * rθφ[1] / (R - outer_radius(s1)) - one(T), 1:legendre_order]

        return [Pr * Yθ for Pr in Prs, Yθ in Ys]
    end

    function kernel_function(r1::T)
        x1 = rθφ2xyz([r1,zero(T),zero(T)])
        Kzero = zeros(Complex{T},len,len_p)

        fun = function (rθφ::SVector{3,T})
            x2 = rθφ2xyz(rθφ)
            if r1 + rθφ[1] < a12
                return Kzero
            end
            if pair_corr(x1,s1,x2,s1) < ptol || norm(x1 - x2) < a12
                return Kzero
            end
            basis2 = field_basis(rθφ)

            U = Uout(x1 - x2)
            U = U .* (bar_numdensity(x2,s1) * pair_corr(x1,s1,x2,s1) * sin(rθφ[2]) * rθφ[1]^2)

            data = [
                sum([
                        t_diags[1][lm_to_n(l,0)] * U[lm_to_n(dl,dm),lm_to_n(l,0)] * b2
                    for b2 in basis2[:,lm_to_n(dl,dm)][:] for l in 0:basis_order]
                for dm = -dl:dl)
            for dl in 0:basis_order];

            return reshape(vcat(data...), (len, len_p))
        end

        return fun
    end

    function δφj(rθφ1::AbstractVector{T})
        basis1 = field_basis(rθφ1)
        return reshape(
            [
                (dl == l) ? b1 : zero(Complex{T})
            for dl in 0:basis_order for b1 in basis1[:,lm_to_n(dl,0)][:] for l in 0:basis_order],
        (len, len_p))
    end

    δφj(r1::T) = δφj([r1,zero(T),zero(T)])

    test_ker = kernel_function(mean(r1s))
    (maxvalue,estimate_error) = hcubature(test_ker, SVector(0.0,0.0,-π), SVector(R-outer_radius(s1),π,π);
        rtol=rtol, maxevals=maxevals
    );

    println("The estimated max coefficient of the integrated kernel is:")
    println(maximum(abs.(maxvalue)) )
    println("with an estimated error of: ")
    println(estimate_error)

    Ks = [
        begin
            ker = kernel_function(r1)
            ker_integrated = hcubature(ker, SVector(0.0,0.0,-π), SVector(R-outer_radius(s1),π,π);
                rtol=rtol, maxevals=maxevals
            )[1]

            δφj(r1) - ker_integrated
        end
    for r1 in r1s];

    bigK = vcat(Ks...);

    bs = incident_coefficients(r1s);
    Fs = bigK \ bs;

    ## Alternative:
    # as = inv(transpose(conj.(bigK)) * bigK) * transpose(conj.(bigK)) * bs;

    if norm(bigK * Fs - bs) / norm(bs) > rtol
        @warn "Numerical solution has a relative residual error of $(norm(bigK * Fs - bs) / norm(bs)), where the requested relative tolernance was: $rtol. This residual error can be decreased by increasing the basis_field_order (current value: $basis_field_order) for the field."
    end

    function scattered_field(x::AbstractVector{T})
        rθφ = cartesian_to_radial_coordinates(x)
        basis1 = field_basis(rθφ)

        Fs = reshape(Fs, (legendre_order,len))

        # F = [[k,dl] for dl in 0:basis_order for k = 1:3]
        # basis1[k,(l,m)] * F[(k,l)]
        return [sum(basis1[:, lm_to_n(l,m)] .* Fs[:,l+1]) for l = 0:basis_order for m = -l:l]
    end

    return ScatteringCoefficientsField(ω, source.medium, material, scattered_field;
        symmetry = RadialSymmetry{3}(),
        basis_order = basis_order,
        basis_field_order = basis_field_order
    )
end

function discrete_system_radial(ω::T, source::AbstractSource{Acoustic{T,3}}, material::Material{3,Sphere{T,3}}, ::RadialSymmetry{3};
        basis_order::Int = 1,
        basis_field_order::Int = Int(round(T(2) * real(ω / source.medium.c) * outer_radius(material.shape))) + 1,
        # basis_field_order = Int(round(T(2) * real(ω / source.medium.c) * outer_radius(material.shape))) + 1,
        polynomial_order::Int = 10,
        legendre_order::Int = 2basis_field_order + 1,
        mesh_points::Int = 3 * (legendre_order + 1),
        # mesh_points = Int(basis_field_order + 1)^2,
        numdensity::Function = (x1, s1) -> number_density(s1),
        # numdensity = (x1, s1) -> number_density(s1),
        scheme = :trapezoidal, pair_kws...
    ) where T

    # basis_field_order = 14

    if length(material.species) > 1
        @warn "discrete_system has only been implemented for 1 species for now. Will use only first specie."
    end

    s1 = material.species[1]
    scale_number_density = one(T) - one(T) / material.numberofparticles
    bar_numdensity(x1,s1) = scale_number_density * numdensity(x1,s1)

    a12 = 2.0 * s1.exclusion_distance * outer_radius(s1)
    R = outer_radius(material.shape)
    k = ω / source.medium.c

    function calculate_gls_fun(;
            pair_corr::Function = hole_correction_pair_correlation,
            pair_corr_distance::Function = z ->  pair_corr(
                [0.0,0.0,0.0],material.species[1], [0.0,0.0,z],material.species[1]
            ),
            sigma_approximation = false,
            gls_pair_radial::Function = gls_pair_radial_fun(
                pair_corr_distance, a12;
                polynomial_order = polynomial_order,
                sigma_approximation = sigma_approximation
            ),
            h12 = R
        )
        gls_pair_radial, h12
    end

    gls_function, h12 = calculate_gls_fun(; pair_kws...)

    ls = 0:basis_order
    lm_to_n = lm_to_spherical_harmonic_index
    ls_to_ns = lm_to_n.(ls,0)

    t_matrices = get_t_matrices(source.medium, material.species, ω, basis_order)
    t_diags = diag.(t_matrices)

    rs = LinRange(0, R - outer_radius(s1), mesh_points)
    σs = integration_scheme(rs; scheme = scheme)

    # incident wave coefficients
    b0 = regular_spherical_coefficients(source)(1,origin(material.shape),ω)[1];
    # regular_radial_basis(source.medium, ω, basis_order, r1)[ls_to_ns]

    Bs = (sqrt(4pi) * b0) .* [
        t_diags[1][lm_to_n(l,0)] * (-T(1))^l * sbesselj(l, k*r1)
    for l in ls, r1 in rs][:]

    chi(l5::Int,l6::Int,r1::T,r2::T) = T(-1)^l6 *
        if(r1 < r2)
            shankelh1(l6, k*r2) * sbesselj(l5, k*r1)
        else
            shankelh1(l5, k*r1) * sbesselj(l6, k*r2)
        end

    function C_kernal(l2::Int,j2::Int)

        term2 = σs[j2] * rs[j2]^2 * bar_numdensity([T(0),T(0),rs[j2]],s1)
        a12 = 2*outer_radius(s1)*s1.exclusion_distance

        data = term2 .* [
        begin
            t_diags[1][lm_to_n(l,0)] *
            if r1 + rs[j2] < a12
                zero(Complex{T})
            elseif abs(r1 - rs[j2]) > h12
                chi(l,l2,r1,rs[j2]) * 4pi * T(-1)^l * T(2l2 + 1)
            else
                gls = gls_function(r1,rs[j2])
                L1 = polynomial_order
                sum(
                    gls[l1+1] *
                    sum(
                        (chi(l5,l6,r1,rs[j2]) * Complex{T}(im)^(l-l5-l6-l2) /
                        (4pi)^2) *
                        sum(
                            T(-1)^m2 *
                            gaunt_coefficient(l,0,l1,m1,l5,-m1) * gaunt_coefficient(l1,m1,l2,m2,l6,m1-m2) *
                            gaunt_coefficient(l2,m2,l,0,l4,m2) * gaunt_coefficient(l4,m2,l5,m1,l6,m2-m1)
                        for m1 = max(-l1,-l5):min(l1,l5)
                        for m2 = max(-l2,-l4,m1-l6):min(l2,l4,m1+l6))
                    for l4 = abs(l-l2):(l+l2)
                    for l5 = abs(l-l1):(l+l1)
                    for l6 = max(abs(l4-l5),abs(l1-l2)):min(l4+l5,l1+l2))
                for l1 = 0:L1)
            end
        end
        for l in ls, r1 in rs];

        return data[:]
    end

    Cs = [C_kernal(l2,j2) for l2 in ls, j2 in eachindex(rs)][:];

    bigC = hcat(Cs...);
    #NOTE bigC[:,M2] == Cs[M2]

    A = I - bigC

    Fs = A \ Bs;

    # NOTE: G * Fs = 0

    # CHECK:
    # norm(Fs - (Bs + bigC * Fs)) / norm(Bs)
    # norm(Fs - (Bs + bigC * Fs))

    # Fc = reshape(Fc,(basis_order+1,length(rs)))
    Fs = reshape(Fs,(basis_order+1,length(rs)))
    # reshape([(l,r) for l in ls, r in rs][:],(basis_order+1,length(rs)))

    # Approximate with a Legendre series
    P = Legendre{T}()

    r1_max = maximum(rs);
    rbars = T(2.0) .* rs ./ r1_max .- T(1.0);
    Pmat = P[rbars, 1:(legendre_order + 1)];

    # Pmat * pls_arr ~ tranpose(Fs)
    pls_arr = transpose(Pmat \ transpose(Fs))

    FsP = pls_arr * transpose(Pmat)
    # Fs ~ pls_arr * transpose(Pmat)

    if  norm(A * FsP[:] - Bs)/norm(Bs) > 0.001
        @warn "Numerical solution has a relative residual error of $(norm((I - bigC) * FsP[:] - Bs)/norm(Bs)), which is larger than expected. This residual error can be decreased by increasing the basis_field_order (current value: $basis_field_order) for the field."
    end

    function scattered_field(xs::AbstractVector{T})
        r,θ,φ = cartesian_to_radial_coordinates(xs)

        Ps = P[2r / r1_max - T(1), 1:(legendre_order + 1)]
        Fls = pls_arr * Ps

        # Fls_arr = [
        #     pls_arr * P[2r / r1_max - T(1), 1:(legendre_order + 1)]
        # for r in rs]
        # norm(Fs - hcat(Fls_arr...)) / norm(Fs)

        conj_Ys = conj.(spherical_harmonics(basis_order, θ, φ));
        fns = [
            Fls[l+1] * conj_Ys[lm_to_n(l,m)]
        for l = 0:basis_order for m = -l:l]

        return fns
    end

    # NOTE: need to think about this type, it is exactly the same as the discrete method above
    return ScatteringCoefficientsField(ω, source.medium, material, scattered_field;
        symmetry = RadialSymmetry{3}(),
        basis_order = basis_order,
        basis_field_order = basis_field_order
    )
end


"""
    outgoing_translation_matrix(ω, ::Acoustic, material::Material{Dim,Sphere{T,Dim}};
        basis_order = 2, tol = 1e-3)

    return a function U where U(X) for X ∈ material.shape gives the the outgoing translation matrix
"""
function outgoing_translation_matrix(ω::T, medium::Acoustic{T,Dim}, material::Material{Dim,Sphere{T,Dim}};
        basis_order::Int = 2, tol::T = 1e-3) where {T,Dim}

    L = basisorder_to_basislength(typeof(medium), basis_order);
    R = outer_radius(material.shape) - minimum(outer_radius.(material.species))
    a12 = T(2) * minimum(outer_radius(s) * s.exclusion_distance for s in material.species)

    rθφ2xyz = radial_to_cartesian_coordinates

    # Just to get the type Tinter.
        xs = LinRange(a12,2R,4);
        data = [zeros(Complex{T},L,L) for x in xs];
        inter = interpolate(
            [data[j][1] for j in CartesianIndices(data)],
            BSpline(Cubic(Line(OnGrid())))
        );
        Tinterx = typeof(scale(inter, xs))

        data = [ zeros(Complex{T},L,L) for x in xs, y in xs, z in xs];
        inter = interpolate(
            [data[j][1] for j in CartesianIndices(data)],
            BSpline(Cubic(Line(OnGrid())))
        );
        Tinterxyz = typeof(scale(inter, xs, xs, xs))

    # Determine the number of mesh points
        Ns = 10:4:66;
        Nθs = 10:4:62;

        xs = LinRange(a12,2R,Ns[end]);
        θs = LinRange(zero(T),T(pi),Nθs[end]);

        datax = [
            outgoing_translation_matrix(medium, basis_order, ω, [x,zero(T),zero(T)])
        for x in xs];

        dataθ = [
            outgoing_translation_matrix(medium, basis_order, ω, rθφ2xyz(SVector(a12,θ,zero(T))))
        for θ in θs];

        CI1s = CartesianIndices(datax[1]);
        inter_xs = Array{Tinterx}(undef, size(CI1s));
        inter_θs = Array{Tinterx}(undef, size(CI1s));

        for i in CI1s
            inter = interpolate(
                [datax[j][i] for j in CartesianIndices(datax)],
                BSpline(Cubic(Line(OnGrid())))
            );
            inter_xs[i] = scale(inter, xs)

            inter = interpolate(
                [dataθ[j][i] for j in CartesianIndices(dataθ)],
                BSpline(Cubic(Line(OnGrid())))
            );
            inter_θs[i] = scale(inter, θs)
        end

        errorsx = Array{T}(undef, length(Ns))
        for n in eachindex(Ns)
            xNs = LinRange(a12,2R,Ns[n]);
            inter_xNs = Array{Tinterx}(undef, size(CI1s));

            for i in CI1s
                inter = interpolate(inter_xs[i].(xNs),BSpline(Cubic(Line(OnGrid()))));
                inter_xNs[i] = scale(inter, xNs)
            end

            errorsx[n] = norm(norm(inter_xNs[i].(xs) - inter_xs[i].(xs)) for i in CI1s) / norm(norm(inter_xs[i].(xs)) for i in CI1s)
        end
        n = findfirst(errorsx .< tol)
        Nr = Ns[n]

        errorsθ = Array{T}(undef, length(Nθs))
        for n in eachindex(Nθs)
            θNs = LinRange(zero(T),T(pi),Nθs[n]);
            inter_θNs = Array{Tinterx}(undef, size(CI1s));

            for i in CI1s
                inter = interpolate(inter_θs[i].(θNs),BSpline(Cubic(Line(OnGrid()))));
                inter_θNs[i] = scale(inter, θNs)
            end

            errorsθ[n] = norm(norm(inter_θNs[i].(θs) - inter_θs[i].(θs)) for i in CI1s) / norm(norm(inter_θs[i].(θs)) for i in CI1s)
        end
        n = findfirst(errorsθ .< tol)
        Nθ = Nθs[n]

    # interpolate the whole spherical region based on the xs mesh size
    rs = LinRange(a12,2R,Nr);
    φs = LinRange(-T(pi),T(pi),Nθ);
    θs = LinRange(zero(T),T(pi),Nθ);

    data = [
            outgoing_translation_matrix(medium, basis_order, ω, rθφ2xyz(SVector(r,θ,φ)))
    for r in rs, θ in θs, φ in φs];

    # reorganise the data to interpolate each element of outgoing_translation_matrix in terms of x, y, z.
    CI1s = CartesianIndices(data[1]);
    inters = Array{Tinterxyz}(undef, size(CI1s));

    for i in CI1s
        inter = interpolate(
            [data[j][i] for j in CartesianIndices(data)],
            BSpline(Cubic(Line(OnGrid())))
        );
        inters[i] = scale(inter, rs, θs, φs)
    end

    function Uout(X::AbstractVector{T})
        if norm(X) < a12
            zeros(Complex{T},L,L)
        else
            rθφ = cartesian_to_radial_coordinates(X)
            map(CI1s) do i
                inters[i](rθφ[1], rθφ[2], rθφ[3])
            end
        end
    end

    return Uout
end
