# Here we calculate the effective wavenumber and effective wave amplitudes, without any restriction to the volume fraction of particles and incident wave frequency.

function effective_wavenumber{T<:Number}(ω::T, medium::Medium{T}, species::Array{Specie{T}}; hankel_order::Int = 3, radius_multiplier = 1.005, kws...)
    k = ω/medium.c
    @memoize Z_l_n(l,n) = Zn(ω,species[l],medium,n)
    @memoize Nn(n,x,y) = x*diffhankelh1(n,x)*besselj(n,y) - y*hankelh1(n,x)*diffbesselj(n,y)

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    M(kef,j,l,m,n) = (n==m ? 1:0)*(j==l ? 1:0) + 2.0pi*species[l].num_density*Z_l_n(l,n)*
            Nn(n-m,k*as[j,l],kef*as[j,l])/(k^2.0-kef^2.0)
    ho = hankel_order
    S = length(species)
    # this matrix is needed to calculate the eigenvectors
    MM(kef) = reshape(
        [ M(kef,j,l,m,n) for m in -ho:ho, j = 1:S, n in -ho:ho, l = 1:S]
    , ((2ho+1)*S, (2ho+1)*S))
    # det(reshape(MM(kef), ( (2ho+1)*length(species), (2ho+1)*length(species)) ))
end
# Zn{T}(ω, p::Specie{T}, med::Medium{T},  m::Int)

ωs = collect(linspace(0.01,1.0,60))
ω = 0.9
species = [
    Specie(ρ=WaterDistilled.ρ,r=30.e-6, c=WaterDistilled.c, volfrac=0.1),
    Specie(ρ=Inf, r=100.0e-6, c=2.0, volfrac=0.2)
]
# background medium
medium = Glycerol
hankel_order = 3
radius_multiplier = 1.005

effective_wavenumber(ω, medium, species)
