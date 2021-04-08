using EffectiveWaves, LinearAlgebra
using HCubature
using Interpolations
import StaticArrays: SVector
import MultipleScattering: outgoing_translation_matrix

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

    carts = radial_to_cartesian_coordinates

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
            outgoing_translation_matrix(medium, basis_order, ω, carts(SVector(a12,θ,zero(T))))
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
            outgoing_translation_matrix(psource.medium, basis_order, ω, carts(SVector(r,θ,φ)))
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


"""
    azimuthal_integral(k,r1_vec::Vector,r2,θ2; basis_order::Int = 3)

documentation
"""
function discrete_system(ω::T, psource::PlaneSource{T,Dim,1,Acoustic{T,Dim}}, material::Material{Dim,Sphere{T,Dim}};
        basis_order::Int = 2,
        rtol = 1e-4, atol=1e-5, maxevals=Int(1e4),
        pair_corr = hole_correction_pair_correlation
    ) where {T,Dim}

    rads = cartesian_to_radial_coordinates
    carts = radial_to_cartesian_coordinates

    k = 0.4;
    ω = 0.4;

    T = Float64

    basis_order = 1
    Dim = 3
    medium = Acoustic(Dim; ρ=1.0, c=1.0)
    psource = PlaneSource(medium, [0.0,0.0,1.0]);

    s1 = Specie(
       Acoustic(3; ρ=10.2, c=10.1), 0.5;
       volume_fraction=0.2
    );

    sphere = Sphere(5.0)
    material = Material(sphere,[s1])

    ω = 0.4;

    a12 = 2.0;
    r1s = 0.1:0.1:1.0;
    θ1s = 0.1:0.2:π;

    R = 5.0
    basis_order = 1;

    rtol = 1e-4; atol = 1e-5; maxevals=Int(1e4);

    Uout = outgoing_translation_matrix(ω, psource.medium, material; basis_order = basis_order, tol = rtol * 10)

    function kernel_function(p::Int, rθ1::SVector{2,T})
        x1 = carts(SVector(rθ1[1],rθ1[2],zero(T)))

        fun = function (rθ2::SVector{3,T})
            x2 = carts(rθ2)
            if norm(x1 - x2) <= a12
                L = basisorder_to_basislength(Acoustic{T,Dim}, basis_order)
                return zeros(Complex{T},L,L)
            else
                # U = outgoing_translation_matrix(psource.medium, basis_order, ω, x1 - x2)
                return Uout(x1 - x2) .* (sin(rθ2[3]) * rθ2[1]^2)
            end
        end

        return fun
    end

    r1 = 3.0; θ1 = 1.3;
    x1 = carts(SVector(r1,θ1,zero(T)))
    ker = kernel_function(1,SVector(r1,θ1))

    rθφ2 = SVector(0.1,0.3,-1.30);
    x2 = carts(rθφ2)
    ker(rθφ2)


    r1 = 3.0; θ1 = 1.3;
    x1 = carts(SVector(r1,θ1,zero(T)))
    # rtol = 1e-4; atol=1e-5; maxevals=Int(1e4);

    # @time outgoing_translation_matrix(psource.medium, basis_order, ω, rand(3))


    r1s = 0.1:0.1:1.0;
    θ1s = 0.1:0.2:π;


    data = [
        begin
            ker = kernel_function(rand(1:5),SVector(r1,θ1))
            hcubature(ker, SVector(0.0,0.0,-pi), SVector(R-1.0,π,π); rtol=rtol, atol=atol, maxevals=maxevals)
        end
    for r1 in r1s, θ1 in θ1s]

end

#
# @time discrete_system(ω, psource, material; basis_order = 1)
# @time discrete_system(ω, psource, material; basis_order = 2)
# # @time discrete_system(ω, psource, material; basis_order = 3)
