"Calculate the largest needed order for the hankel series."
function maximum_hankel_order(ω::Union{T,Complex{T}}, medium::Medium{T}, species::Vector{Specie{T}};
        tol::T=1e-7, verbose::Bool = false) where T <: Number

    # estimation is based on far-field scattering pattern contribution to the primary effective wavenumber^2 and reflection coefficient.
    f0 = ω^(2.0)/(medium.c^2)
    hankel_order = 0
    next_order = -4im*sum(sp.num_density*Zn(ω,sp,medium,hankel_order) for sp in species)

    hankel_order = 1
    next_order = -4im*sum(sp.num_density*Zn(ω,sp,medium,hankel_order) for sp in species)

    # increase hankel order until the relative error^2 < tol. Multiplying the tol by 200 has been chosen based on empircal comparisons with other tolerance used for reflection coefficients.
    while abs(next_order/f0) > tol*200
        f0 = f0 + next_order
        hankel_order += 1
        next_order = -4im*sum(sp.num_density*Zn(ω,sp,medium,hankel_order) for sp in species)
    end

    return hankel_order
end

"A t_matrix in the form of a vector, because for now we only deal with diagonal T matrices."
function t_vectors(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; hankel_order = 3, dim = 2) where T <: AbstractFloat
    t_vecs = [ zeros(Complex{T},1+2hankel_order) for s in species]
    for i = 1:length(species), n = 0:hankel_order
        t_vecs[i][n+hankel_order+1] = - Zn(ω,species[i],medium,n; dim = dim)
        t_vecs[i][-n+hankel_order+1] = t_vecs[i][n+hankel_order+1]
    end
    return t_vecs
end

"Pre-calculate a matrix of Zn's"
function Zn_matrix(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; hankel_order = 3) where T <: Number
    Zs = OffsetArray{Complex{T}}(undef, 1:length(species), -hankel_order:hankel_order)
    for i = 1:length(species), n = 0:hankel_order
        Zs[i,n] = Zn(ω,species[i],medium,n)
        Zs[i,-n] = Zs[i,n]
    end
    return Zs
end

Zn(ω::T, p::Specie{T}, med::Medium{T}, m::Int; kws...) where T<: Number = Zn(Complex{T}(ω), p, med, m; kws...)

"Returns a ratio which gives the scattering strength of a particle"
function Zn(ω::Complex{T}, p::Specie{T}, med::Medium{T},  m::Int; dim::Int = 2) where T<: Number
    m = abs(m)
    ak = p.r*ω/med.c
    # check for material properties that don't make sense or haven't been implemented
    if abs(p.c*p.ρ) == T(NaN)
        error("scattering from a particle with density =$(p.ρ) and phase speed =$(p.c) is not defined")
    elseif abs(med.c*med.ρ) == T(NaN)
        error("wave propagation in a medium with density =$(med.ρ) and phase speed =$(med.c) is not defined")
    elseif abs(med.c) == zero(T)
        error("wave propagation in a medium with phase speed =$(med.c) is not defined")
    elseif abs(med.ρ) == zero(T) && abs(p.c*p.ρ) == zero(T)
        error("scattering in a medium with density $(med.ρ) and a particle with density =$(p.ρ) and phase speed =$(p.c) is not defined")
    end
    if dim == 3
        j = sbesselj
        dj(m,x) = diffsbessel(sbesselj,m,x)
        h1 = shankelh1
        dh1(m,x) = diffsbessel(shankelh1,m,x)
    else
        j = besselj
        dj = diffbesselj
        h1 = hankelh1
        dh1 = diffhankelh1
    end
    # set the scattering strength and type
    if abs(p.c) == T(Inf) || abs(p.ρ) == T(Inf) || abs(med.ρ) == zero(T)
        numer = dj(m, ak)
        denom = dh1(m, ak)
    else
        q = (p.c*p.ρ)/(med.c*med.ρ) #the impedance
        if q == zero(T)
          numer =  j(m, ak)
          denom =  h1(m, ak)
        else
          γ = med.c/p.c #speed ratio
          numer = q * dj(m, ak)*j(m, γ*ak) - j(m, ak)*dj(m, γ*ak)
          denom = q * dh1(m, ak)*j(m, γ*ak) - h1(m, ak)*dj(m, γ*ak)
        end
    end

    return numer / denom
end
