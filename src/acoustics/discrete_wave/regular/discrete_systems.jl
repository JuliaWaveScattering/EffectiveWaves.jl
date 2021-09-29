using EffectiveWaves, LinearAlgebra
using HCubature
using Interpolations, ClassicalOrthogonalPolynomials
import StaticArrays: SVector

import MultipleScattering: outgoing_translation_matrix

"""
    discrete_system(k,r1_vec::Vector,r2,θ2; basis_order::Int = 3)

documentation
"""
function discrete_system(ω::T, source::Source{T,Acoustic{T,Dim}}, material::Material{Dim,Sphere{T,Dim}}, ::AbstractAzimuthalSymmetry{Dim};
        basis_order::Int = 2,
        field_basis_order::T = 4
        rtol = 1e-3, atol=1e-3, maxevals=Int(2e4),
        pair_corr = hole_correction_pair_correlation
    ) where {T,Dim}

    pair_corr = hole_correction_pair_correlation

    ω = 0.01;
    T = Float64

    rtol = 1e-3; atol = 1e-3; maxevals=Int(2e4);

    field_basis_order = 4;
    basis_order = 1
    Dim = 3
    medium = Acoustic(Dim; ρ=1.0, c=1.0)
    psource = PlaneSource(medium, [0.0,0.0,1.0]);
    source = plane_source(medium; direction = [0.0,0.0,1.0])

    s1 = Specie(
       Acoustic(3; ρ=10.2, c=10.1), 0.5;
       volume_fraction=0.2
    );

    sphere = Sphere(5.0)
    material = Material(sphere,[s1])

    scale_number_density = one(T) - one(T) / material.numberofparticles
    bar_numdensity = scale_number_density * number_density(s1)

    R = outer_radius(material.shape)

    r1s = LinRange(0,R-outer_radius(s1),Int(round(1.1*field_basis_order)) + 1)
    θ1s = LinRange(0,π,Int(round(1.1*field_basis_order)) + 1)

    Uout = outgoing_translation_matrix(ω, psource.medium, material; basis_order = basis_order, tol = atol * 10)

    P = Legendre{T}()

    function field_basis(rθ::AbstractVector{T})
        # [P_0(cos(θ)), …, P_(field_basis_order-1)(cos(θ))]
        Pθs = P[cos(rθ[2]),1:field_basis_order]
        Prs = P[2 * rθ[1] / R - one(T),1:field_basis_order]

        # [Pr * Pθ for Pr in Prs, Pθ in Pθs][:]
        return (Prs * transpose(Pθs))[:]
    end

    L = basisorder_to_basislength(Acoustic{T,Dim}, basis_order)
    ls, ms = spherical_harmonics_indices(basis_order)

    a12s = [
        outer_radius(s1) * s1.exclusion_distance + outer_radius(s2) * s2.exclusion_distance
    for s1 in material.species, s2 in material.species]

    t_matrices = get_t_matrices(psource.medium, material.species, ω, basis_order)
    t_diags = diag.(t_matrices)

    rθφ2xyz = radial_to_cartesian_coordinates

    function kernel_function(rθ1::SVector{2,T})
        x1 = rθφ2xyz(SVector(rθ1[1],rθ1[2],zero(T)))

        fun = function (rθφ::SVector{3,T})
            x2 = rθφ2xyz(rθφ)
            if norm(x1 - x2) <= a12s[1,1]
                return zeros(Complex{T},L,L*field_basis_order^2)
            end
            basis1 = field_basis(SVector(rθφ[1],rθφ[2]))
            U = Uout(x1 - x2)
            U = U .* (bar_numdensity * pair_corr(x1,s1,x2,s1) * sin(rθφ[2]) * rθφ[1]^2)
            return reshape(
                [
                    ((nd == n) ? b : zero(Complex{T})) - t_diags[1][n] * U[nd,n] * b * exp(-im*ms[nd]*rθφ[3])
                for n in 1:L, nd in 1:L, b in basis1],
            (L,L*field_basis_order^2))
        end

        return fun
    end

    data = [
        begin
            ker = kernel_function(SVector(r1,θ1))
            hcubature(ker, SVector(0.0,0.0,-π), SVector(R-outer_radius(s1),π,π);
                rtol=rtol, atol=atol, maxevals=4*maxevals
                rtol=rtol, atol=atol, maxevals=maxevals
            )[1]
        end
    for r1 in r1s, θ1 in θ1s];

    K = vcat(data...);

    gs = regular_spherical_coefficients(source)(basis_order,origin(material.shape),ω);

    v = regular_basis_function(medium,  ω)
    lm2n = lm_to_spherical_harmonic_index

    bs = [
        begin
            vs = v(field_basis_order, rθφ2xyz(SVector(r1,θ1,zero(T))))
            [
                sum(
                    (abs(m) > l1) ? zero(Complex{T}) :
                        gaunt_coefficient(dl,0,l,m,l1,m) * vs[lm2n(l1,m)] * gs[lm2n(dl,0)]
                for l1 in 0:field_basis_order, dl in 0:basis_order)
            for l = 0:basis_order for m = -l:l]
        end
    for r1 in r1s, θ1 in θ1s];

    bs = vcat(bs...);

    as = inv(transpose(conj.(K)) * K) * transpose(conj.(K)) * bs;
    as = reshape(as,(L,field_basis_order^2));

    fs(rθ::Vector{T}) = as * field_basis(rθ)


    return fs
end

#
# @time discrete_system(ω, psource, material; basis_order = 1)
# @time discrete_system(ω, psource, material; basis_order = 2)
# # @time discrete_system(ω, psource, material; basis_order = 3)


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
        data = [ zeros(Complex{T},L,L) for x in xs];
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
            outgoing_translation_matrix(psource.medium, basis_order, ω, rθφ2xyz(SVector(r,θ,φ)))
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
