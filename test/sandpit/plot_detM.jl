using EffectiveWaves, Plots

pyplot()

Maxtime=100.
T=Float64

"Derivative of Hankel function of the first kind"
function diffhankelh1(n,z)
  if n!=0
    0.5*(hankelh1(-1 + n, z) - hankelh1(1 + n, z))
  else
    - hankelh1(1, z)
  end
end

"Derivative of Bessel function of first kind"
function diffbesselj(n,z)
  if n!=0
    0.5*(besselj(-1 + n, z) - besselj(1 + n, z))
  else
    - besselj(1, z)
  end
end

# background medium
medium = Medium(1.0,1.0+0.0im)
radius_multiplier = 1.005
max_basis_order=15

incident_medium = medium
θin = 0.2
ω = 15.
k = ω/incident_medium.c

## Strong scatterers
species = [
    Specie(ρ=0.8,r=0.1, c=0.2, volfrac=0.1),
    Specie(ρ=0.2, r=0.2, c=0.1, volfrac=0.1)
]

S = length(species)
# @memoize Z_l_n(l,n) = Zn(ω,species[l],medium,n)

as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
function M(keff,j,l,m,n)
    (n==m ? 1.0+im*0.0:0.0+im*0.0)*(j==l ? 1.0+im*0.0:0.0+im*0.0) + 2.0pi*species[l].num_density*Z_l_n[l,n]*
        kernelN(n-m,k*as[j,l],keff*as[j,l])/(k^2.0-keff^2.0)
end

ho = -1 + sum([ tol .< norm([M(0.9*k + 0.1im,j,l,1,n) for j = 1:S, l = 1:S]) for n=0:max_basis_order ])

# this matrix is needed to calculate the eigenvectors
MM(keff::Complex{T}) = reshape(
    [M(keff,j,l,m,n) for m in -ho:ho, j = 1:S, n in -ho:ho, l = 1:S]
, ((2ho+1)*S, (2ho+1)*S))
detMM2(keff_vec::Array{T}) = map(x -> real(x*conj(x)), det(MM(keff_vec[1]+im*keff_vec[2])))

function detMM!(F,x)
    F[1] = abs(det(MM(x[1]+im*x[2])))
end

k0 = real(k)
x = k0 .* LinRange(0.0,1.8,100)/10.
y = k0 .* LinRange(0.,0.22,100)

X = repmat(x',length(y),1)
Y = repmat(y,1,length(x))

Z0 = map( (x,y) -> detMM2([x,y]),X,Y)
Z = map( z -> (abs(z)> 2.) ? NaN : z, Z0)
contour(x,y,Z,fill=true, xlab = "Re k*", ylab = "Im k*", title="Roots of secular det M - high ω - low φ")

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
