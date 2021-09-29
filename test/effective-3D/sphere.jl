using EffectiveWaves
using LinearAlgebra

include("test_discrete_solution.jl")
# Set parameters
particle_medium = Acoustic(3; ρ=0.01, c=0.01);
particle_medium = Acoustic(3; ρ=10.0, c=10.0);
medium = Acoustic(3; ρ=1.0, c=1.0);

R = 8.0
r = 1.0

separation_ratio = 1.02

ka = 0.4
k = ka / r

vol_fraction = 0.1

basis_order = Int(round(4. * ka)) .+ 1
basis_field_order = Int(round(2.0 * k * R)) .+ 1

# basis_order = 1;
# basis_field_order = 6;
# basis_field_order = 10;

ω = k * real(medium.c)

s1 = Specie(
    particle_medium, Sphere(r);
    number_density = vol_fraction / volume(Sphere(r)),
    exclusion_distance = separation_ratio
);
species = [s1]

tol = 1e-7

# eff_medium = effective_medium(medium, species)

psource = PlaneSource(medium, [0.0,0.0,1.0]);
source = plane_source(medium; direction = [0.0,0.0,1.0])

region_shape = Sphere([0.0,0.0,0.0], R)
material = Material(Sphere(R),species);

eff_medium = effective_medium(medium, species; numberofparticles = material.numberofparticles)
k_low = ω / eff_medium.c

keffs = wavenumbers(ω, medium, species;
    # num_wavenumbers = 2,
    mesh_points = 15,
    mesh_size = 2.0,
    num_wavenumbers = 3,
    basis_order = basis_order,
    tol = tol,
    symmetry = PlanarAzimuthalSymmetry(),
    numberofparticles = material.numberofparticles
)

pwavemode = WaveMode(ω, k_low, psource, material;
    basis_order = basis_order,
    basis_field_order = basis_field_order
);

basis_field_order = 6

pwavemode = WaveMode(ω, keffs[1], psource, material;
    basis_order = basis_order,
    basis_field_order = basis_field_order
);

wavemode = WaveMode(ω, keffs[1], source, material;
    basis_order = basis_order,
    basis_field_order = basis_field_order
);

scat_field = scattering_field(wavemode)
pscat_field = scattering_field(pwavemode)

res = discrete_system_residue(pscat_field, ω, source, material, AzimuthalSymmetry{3}();
    basis_order = basis_order, mesh_points = 5,
    rtol = 1e-2, maxevals = Int(1e4)
)

# xs = [ [cos(θ), sin(θ), 1.1] for θ in 0.0:0.6:(2π), r in LinRange(0.0,R,)]
xs = [ radial_to_cartesian_coordinates([R - 4*outer_radius(s1),θ,0.0]) for θ in 0.0:0.1:(π)];
# xs = [xs; [ radial_to_cartesian_coordinates([2.1*r,θ,0.0]) for θ in 0.0:0.1:(π)]];
# xs = [ radial_to_cartesian_coordinates([r,0.0,0.0]) for θ in 0.0:0.1:(R-r)];

pdata = pscat_field.(xs);
data = scat_field.(xs);
maximum(abs.(d[1]) for d in pdata - data)

pdata0 = [d[1] for d in pdata];
pdata1 = [d[3] for d in pdata];

pdata11 = [d[4] for d in pdata];
pdata1m1 = [d[2] for d in pdata];

data0 = [d[1] for d in data];
data1 = [d[3] for d in data];


using Plots
using Statistics


rtol = 1e-2; maxevals = Int(5e4);
legendre_order = basis_field_order+2;

pair_corr = hole_correction_pair_correlation;

discrete_scat = discrete_system(ω, psource, material;
    basis_order = basis_order,
    legendre_order = legendre_order,
    rtol = rtol, maxevals = maxevals
);

scatdata = discrete_scat.(xs);

maximum(norm.(scatdata - pdata) ./ norm.(pdata))
mean(norm.(scatdata - pdata) ./ norm.(pdata))

inds = findall(abs.(mat_coefs_pwave) .> 1e-6)

mat_coefs_pwave = material_scattering_coefficients(pwavemode)
# material_scattering_coefficients(wavemode)

mat_coefs_pdisc = material_scattering_coefficients_discrete(pscat_field, ω, source, material;
    basis_order = basis_order,
    basis_field_order = basis_field_order,
    rtol = 1e-2,
    maxevals = Int(1e4)
)

mat_coefs_disc2 = material_scattering_coefficients_discrete(discrete_scat, ω, source, material;
    basis_order = basis_order,
    basis_field_order = basis_field_order,
    rtol = 1e-2,
    maxevals = Int(1e5)
);

L = length(mat_coefs_pdisc)
norm(mat_coefs_pwave[1:L] - mat_coefs_pdisc) / norm(mat_coefs_pwave[1:L])
norm(mat_coefs_pwave[1:L] - mat_coefs_disc) / norm(mat_coefs_pwave[1:L])
norm(mat_coefs_pwave[1:L] - mat_coefs_disc2) / norm(mat_coefs_pwave[1:L])

[mat_coefs_pwave[inds] mat_coefs_disc2[inds] mat_coefs_disc[inds]]

# Calculate low frequency scattering
    Linc = basis_order + basis_field_order
    source_coefficients = regular_spherical_coefficients(source)(Linc,zeros(3),ω);

    material_low = Material(
        Sphere(outer_radius(material.shape) - outer_radius(s1)),
        species
    );

    effective_sphere = Particle(eff_medium, material_low.shape);
    Tmat = MultipleScattering.t_matrix(effective_sphere, medium, ω, Linc);
    scat_coef_low = Tmat * source_coefficients;


# 0.747
scatdata0 = [d[1] for d in scatdata];
scatdata1 = [d[3] for d in scatdata];
scatdata11 = [d[4] for d in scatdata];
scatdata1m1 = [d[2] for d in scatdata];

plot([real.(scatdata0),real.(pdata0)])
plot([imag.(scatdata0),imag.(pdata0)])

plot([real.(scatdata1),real.(pdata1)])
plot([imag.(scatdata1),imag.(pdata1)])

# discrete_scat = discrete_system(ω, psource, material;
#     basis_order = basis_order, basis_field_order = basis_field_order,
#     rtol = 1e-3, atol=1e-3, maxevals=Int(10e4)
# );
# d2 = discrete_scat([1.0,2.0,1.1])
#
# discrete_scat = discrete_system(ω, psource, material;
#     basis_order = basis_order, basis_field_order = basis_field_order,
#     rtol = 1e-4, atol=1e-4, maxevals=Int(20e4)
# );
# d3 = discrete_scat([1.0,2.0,1.1])
# [d1,d2,d3]

plot!(fun.(scatdata0 ))

fun = imag
fun = real
plot([fun.(scatdata0),fun.(pdata0)])


pdata[1]
data[1]

# Test chosen basis for field

    R = outer_radius(material.shape)
    legendre_order = basis_field_order
    mesh_points = Int(round(1.5*legendre_order)) + 1

    r1s = LinRange(0,R-r, mesh_points)
    θ1s = LinRange(0,π, mesh_points)

    T = Float64
    using ClassicalOrthogonalPolynomials
    function field_basis(rθ::AbstractVector{T},legendre_order::Int)
        P = Legendre{T}()

        # [P_0(cos(θ)), …, P_(legendre_order-1)(cos(θ))]
        Pθs = P[cos(rθ[2]), 1:legendre_order]
        Prs = P[2 * rθ[1] / (R - r) - one(T), 1:legendre_order]

        # [Pr * Pθ for Pr in Prs, Pθ in Pθs][:]
        return (Prs * transpose(Pθs))[:]
    end

    field_basis([1.0,0.2])

    len = basisorder_to_basislength(Acoustic{T,3}, basis_order)
    len_p = legendre_order^2
    #
    # function δφj(rθ1::AbstractVector{T}, legendre_order = 5)
    #     basis1 = field_basis(rθ1, legendre_order)
    #     return [b1 for b1 in basis1, n in 1:len][:]
    # end

    function δφj(rθ1::AbstractVector{T}, legendre_order::Int)
        basis1 = field_basis(rθ1, legendre_order)
        return reshape(
            [
                (nd == n) ? b1 : zero(Complex{T})
                # (n,nd,b1)
            # for n in -1:1, nd in 1:3, b1 in 0:2],
            for n in 1:len, nd in 1:len, b1 in basis1],
        (len, len * len_p))
    end

    As = [δφj([r1,θ1], legendre_order) for r1 in r1s, θ1 in θ1s][:];
    A = vcat(As...)

    rθφ2xyz = radial_to_cartesian_coordinates
    b = [pscat_field(rθφ2xyz([r,θ,0.0])) for r in r1s, θ in θ1s][:]
    b = vcat(b...)[:]

    as = A \ b;
    as = reshape(as,(len,:));

    function scattered_field(xs::Vector{T})
        rθφ = cartesian_to_radial_coordinates(xs)
        return as * field_basis(rθφ[1:2],legendre_order)
    end

    rθs = [ [(R-r) * rand(), π * rand(),0.0] for i = 1:2000 ];
    xs = rθφ2xyz.(rθs);

    scatdata = scattered_field.(xs);
    pdata = pscat_field.(xs);

    maximum(norm.(scatdata - pdata) ./ norm.(pdata))
    mean(norm.(scatdata - pdata) ./ norm.(pdata))

# Test reshape used in the discrete_system
    # ps = -1:1
    # ns = 1:3
    # nds = 1:2
    # js = -1:2
    #
    # bs = [
    #         [
    #             (n,j)
    #             # sum((n,j))
    #         for n in ns]
    # for j in js];
    #
    # bs = vcat(bs...);
    #
    # a = [
    #         (nd,p)
    #         # nd + p*1.0im
    # for nd in nds, p in ps][:]
    #
    # function Kmat(j)
    #     reshape(
    #         [
    #             (n,j,nd,p)
    #             # Float64(n)^j + Float64(nd)^p * 1.0im + sum((n,j,nd,p)) + ((n == nd && p == j) ? 10.0 + 2im : 0.0)
    #         for n in ns, nd in nds, p in ps],
    #     (length(ns), length(nds) * length(ps)))
    # end
    #
    # data = [
    #         Kmat(j)
    # for j in js];
    #
    # K = vcat(data...);
    # bs = K * a
    #
    # a2 = inv(transpose(conj.(K)) * K) * transpose(conj.(K)) * bs;
    # a3 = K \ bs;
    # maximum(abs.(a - a2))
    # maximum(abs.(a2 - a3))
    #
    # f(x) = cos(x[1]+x[2])
    # hcubature(f, SVector(0.0,0.0), SVector(30pi,10pi);
    #     rtol=rtol, atol=atol, maxevals=100000
    # )

    function δφ(j)
        return reshape(
            [
                (n,j,nd,p)
            for n in ns, nd in nds, p in ps],
        (length(ns), length(nds) * length(ps)))
    end

    data = [
            δφ(j)
    for j in js];

    Kδφ = vcat(data...);



    legendre_order = 15
    legendre_order = 20
    function field_basis(rθ::AbstractVector{T})
        P = Legendre{T}()
        # [P_0(cos(θ)), …, P_(legendre_order-1)(cos(θ))]
        Pθs = P[cos(rθ[2]), 1:legendre_order]
        Prs = P[2 * rθ[1] / (R - outer_radius(s1)) - one(T), 1:legendre_order]

        # [Pr * Pθ for Pr in Prs, Pθ in Pθs][:]
        return (Prs * transpose(Pθs))[:]
    end


    rs = 0.0:0.1:(R - outer_radius(s1))
    θs = 0.0:0.1:π

    associatedlegendre(m) = ((-1)^m*prod(1:2:(2m-1)))*(UltrasphericalWeight((m+1)/2).*Ultraspherical(m+1/2))
    f(rθ) = sbesselj(3, 0.1 * rθ[1]) * associatedlegendre(2)[cos(rθ[2]),3] * exp(2im) + sbesselj(1,0.1 * rθ[1]) * associatedlegendre(1)[cos(rθ[2]),2] * exp(1im)

    # dot(as,field_basis(rθ)) = f(rθ)

    A = vcat([transpose(field_basis([r,θ])) for r in rs, θ in θs][:]...)
    b = [f([r,θ]) for r in rs, θ in θs][:]

    x = A \ b;

    maximum(abs.(A * x - b)) / mean(abs.(b))

    f1 = [ f([rs[10],θ])  for θ in 0.0:0.03:pi]
    f2 = [ sum(field_basis([rs[10],θ]) .* x) for θ in 0.0:0.03:pi]

    sum(field_basis([0.2,0.3]) .* x)
