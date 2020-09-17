using EffectiveWaves, Plots

pyplot()

Maxtime=100.
T=Float64

spatial_dim = 3
medium = Acoustic(spatial_dim; ρ=1.0, c=1.0)

v = 1.0
s1 = Specie(
    Acoustic(spatial_dim; ρ=10.0, c=10.0), Sphere(1.0);
    volume_fraction = v * 3/4
);
s2 = Specie(
    Acoustic(spatial_dim; ρ=0.1, c=0.1), Sphere(1.0);
    volume_fraction = v * 1/4
);
species = [s1,s2]

k_low = ω ./ effective_medium(medium, species).c


## Swap a specie for the medium

    medium = Acoustic(spatial_dim; ρ=10.0, c=10.0)

    v = 1.0
    s2 = Specie(
        Acoustic(spatial_dim; ρ=0.1, c=0.1), Sphere(1.0);
        volume_fraction = v * 1/4
    );
    species = [s2]

    k_low = ω ./ effective_medium(medium, species).c

ω = 0.04187
basis_order = 1
tol = 1e-7
k = ω/medium.c

opts = Dict(
    :symmetry => PlanarAzimuthalSymmetry(),
    :tol => tol,
    # :num_wavenumbers => 4,
    # :mesh_size => 2.0, :mesh_points => 20,
    :basis_order => basis_order
)

detMM = dispersion_complex(ω, medium, species; opts...)

kps =  wavenumbers(ω, medium,  species; opts...)

k_low = ω ./ effective_medium(medium, species).c

# Explore around all the wavenumbers

inds = 1:1
x_max = maximum(abs.(real.(kps[inds])))
y_max = maximum(abs.(imag.(kps[inds])))

x = LinRange(-x_max,x_max,150)
y = LinRange(0.,y_max,100)

X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(y, 1, length(x))

z_max = 2.0
Z0 = map( (x,y) -> abs(detMM(x + y*im)),X,Y);
Z = map( z -> (abs(z)> z_max) ? NaN : z, Z0);
contour(x,y,Z,fill=true, xlab = "Re k*", ylab = "Im k*", title="Roots of secular det M - high ω - low φ")
scatter!(kps[inds])


# Explore around one wavenumber

    k1 = abs(real(kps[1])) - imag(kps[1]) * im

    scale = 0.01
    x = LinRange(real(k1) - abs(real(k1)) * scale, real(k1) + abs(real(k1)) * scale,100)
    y = LinRange(imag(k1) - abs(imag(k1)) * scale, imag(k1) + abs(imag(k1)) * scale,100)

    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))

    Z0 = map( (x,y) -> abs(detMM(x + y*im)),X,Y);

    z_max = 2 * scale
    Z = map( z -> (abs(z)> z_max) ? NaN : z, Z0);
    contour(x,y,Z,fill=true, xlab = "Re k*", ylab = "Im k*", title="Roots of secular det M - high ω - low φ")
    scatter!([k1])


# Old code below



x = k0 .* LinRange(0.18,0.23,50)
y = k0 .* LinRange(0.01,0.17,50)

x = kx; y = ky;

X = repmat(x',length(y),1)
Y = repmat(y,1,length(x))
Z = map( (x,y) -> (z = detMM2([x,y]); (abs(z)> 0.5) ? NaN: z),X,Y)

contour(x,y,Z,fill=true, xlab = "Re k*", ylab = "Im k*", title="Roots of secular det M")

x = LinRange(3.,6.,40)
y = LinRange(-1.5,14.,40)

X = repmat(x',length(y),1)
Y = repmat(y,1,length(x))
Z = map( (x,y) -> (z = detMM2([x,y]); (abs(z)> 1.) ? NaN: z),X,Y)
# Z = map( (x,y) -> detMM2([x,y]),X,Y)

# contour(x,y,Z,fill=true, xlab = "Re k*", ylab = "Im k*", title="Roots of secular det M")
heatmap(x,y,Z, xlab = "Re k*", ylab = "Im k*", title="det M - high ω and almost dirichlet")
