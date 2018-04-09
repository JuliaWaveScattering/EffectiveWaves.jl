# Here we calculate the effective wavenumber and effective wave amplitudes, without any restriction to the volume fraction of particles and incident wave frequency.

Nn(n,x,y) = x*diffhankelh1(n,x)*besselj(n,y) - y*hankelh1(n,x)*diffbesselj(n,y)

function reflection_coefficient{T<:Number}(ω::T, medium::Medium{T}, species::Array{Specie{T}}; θ_inc::T = 0.0, kws...)
    (k_eff,θ_eff,As) = effective_planewave(ω, medium, species; θ_inc=θ_inc, kws...)
    θ_ref = pi - θ_eff - θ_inc
    k = ω/medium.c
    S = length(species); ho = Int((size(As,1)-1)/2)

    R = 2.0im/(k*cos(θ_inc)*(k*cos(θ_inc) + k_eff*cos(θ_eff)))
    R = R*sum(
            exp(im*n*θ_ref)*species[l].num_density*As[n+ho+1,l]*Zn(ω,species[l],medium,n)
    for n=-ho:ho, l=1:S)
    return R
end

function transmitted_planewave{T<:Number}(ω::T, medium::Medium{T}, species::Array{Specie{T}}; hankel_order = :auto, max_hankel_order = 10,
        radius_multiplier = 1.005,
        MaxTime=100., tol = 1e-6, θ_inc::T = 0.0,
        kws...)
    k = ω/medium.c
    S = length(species)
    @memoize Z_l_n(l,n) = Zn(ω,species[l],medium,n)

    as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
    M(keff,j,l,m,n) = (n==m ? 1.0+im*0.0:0.0+im*0.0)*(j==l ? 1.0+im*0.0:0.0+im*0.0) + 2.0pi*species[l].num_density*Z_l_n(l,n)*
            Nn(n-m,k*as[j,l],keff*as[j,l])/(k^2.0-keff^2.0)

    if hankel_order == :auto
        ho = -1 + sum([ tol .< norm([M(0.9*k + 0.1im,j,l,1,n) for j = 1:S, l = 1:S]) for n=0:max_hankel_order ])
    else ho = hankel_order
    end

    # this matrix is needed to calculate the eigenvectors
    MM(k::Complex{T}) = reshape(
        [M(k,j,l,m,n) for m in -ho:ho, j = 1:S, n in -ho:ho, l = 1:S]
    , ((2ho+1)*S, (2ho+1)*S))
    detMM2(k::Array{T}) = map(x -> real(x*conj(x)), det(MM(k[1]+im*k[2])))

    initial_k_eff = multispecies_wavenumber(ω, medium, species);
    initial_k_eff = [real(initial_k_eff), imag(initial_k_eff)]
    lower = [0.,-1.]
    upper = [1.0, 1.0]*100.0*ω/minimum(map(s -> real(s.c), species))
    result = optimize(detMM2, initial_k_eff, lower, upper)

    # Check result
    k_eff = result.minimizer[1] + im*result.minimizer[2]
    MM_svd = svd(MM(k_eff))
    if last(MM_svd[2]) > tol
        warn("Local optimisation was unsucessful at finding an effective wavenumber ( $(last(MM_svd[2])) was the smallest singular value of the effective wavenumber matrix equation ).")
        println("Note that there is not a unique effective wavenumber. BlackBoxOptim (global minimization package) often picks out another root that leads to strange transmission angles and large amplitudes");
        # result = bboptimize(detMM2; MaxTime = min(100.0,MaxTime), SearchRange = [(0., upper[1]), (-1.0, upper[1])], NumDimensions=2, TargetFitness = 1e-10)
        # k_eff = best_candidate(result)[1] + im*best_candidate(result)[2]
    end

    # calculate effective transmission angle
    snell(θ::Array{T}) = norm(k*sin(θ_inc) - k_eff*sin(θ[1] + im*θ[2]))
    result = optimize(snell, [θ_inc,0.], NelderMead())
    θ_eff = result.minimizer[1] + im*result.minimizer[2]

    # calculate effective amplitudes
    A_null = MM_svd[3][:,(2ho+1)*S] # norm(MM(kef)*A_null) ~ 0
    A_null = reshape(A_null, (2*ho+1,S)) # A_null[:,j] = [A[-ho,j],A[-ho+1,j],...,A[ho,j]]

    sumAs = 2*sum(
            exp(im*n*(θ_inc - θ_eff))*Z_l_n(l,n)*species[l].num_density*A_null[n+ho+1,l]
    for n = -ho:ho, l = 1:S)
    x = im*k*cos(θ_inc)*(k_eff*cos(θ_eff) - k*cos(θ_inc))/sumAs

    (k_eff,θ_eff,A_null*x)
end


using EffectiveWaves, Memoize, SpecialFunctions, Optim, BlackBoxOptim

Maxtime=100.
T=Float64
ωs = collect(linspace(0.01,1.0,60))
ω = 0.9
species = [
    Specie(ρ=10.,r=0.1, c=12., volfrac=0.1),
    Specie(ρ=3., r=0.2, c=2.0, volfrac=0.3)
]
# background medium
medium = Medium(1.0,1.0+0.0im)
hankel_order = 3
radius_multiplier = 1.005

# effective_wavenumber(ω, medium, species)
kef0 = map(x -> [abs(real(x)),abs(imag(x))], multispecies_wavenumber(ω, medium, species))


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
