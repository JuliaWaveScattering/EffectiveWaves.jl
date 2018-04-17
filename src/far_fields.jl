# Before far field patterns and the pair field patterns are calculated here. These are needed for the effective wavenumber and reflection for low volume fraction.

d(x,m) = diffbesselj(m,x)*diffhankelh1(m,x) + (1.0 - (m/x)^2)*besselj(m,x)*hankelh1(m,x)

function far_field_pattern(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; tol=1e-6,
        hankel_order = maximum_hankel_order(ω, medium, species; tol=tol),
        verbose = false) where T<:Number

    if verbose
        println("$hankel_order was the largest hankel order used for the far field pattern")
    end
    Zs = Zn_matrix(ω, medium, species; hankel_order = hankel_order)
    num_density_inv = sum(sp.num_density for sp in species)^(-1.0)

    far_field(θ::T) where T <: Number = -num_density_inv*sum(
        species[i].num_density*Zs[i,n]*exp(im*θ*n)
    for i in eachindex(species), n=-hankel_order:hankel_order)

    return far_field
end

function pair_field_pattern(ω::T, medium::Medium{T}, species::Vector{Specie{T}}; tol = T(1e-6),
        hankel_order = maximum_hankel_order(ω, medium, species; tol=tol),
        radius_multiplier = T(1.005),
        verbose = false) where T<:Number

    Zs = Zn_matrix(ω, medium, species; hankel_order = hankel_order)
    num_density_inv = sum(sp.num_density for sp in species)^(-1.0)

    pair_field(θ::T) where T <: Number = -T(π)*num_density_inv^(2.0)*sum(
        begin
            a12 = radius_multiplier*(species[i].r + species[j].r)
            species[i].num_density*species[j].num_density*a12^2.0*d(ω/medium.c*a12,m-n)*
            Zs[i,n]*Zs[j,m]*exp(im*m*θ)
        end
    for i=1:length(species), j=1:length(species),
    n = -hankel_order:hankel_order, m = -hankel_order:hankel_order)

    return pair_field
end
