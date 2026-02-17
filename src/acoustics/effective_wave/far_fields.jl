# Before far field patterns and the pair field patterns are calculated here. These are needed for the effective wavenumber and reflection for low volume fraction.

d2D(x,m) = diffbesselj(m,x)*diffhankelh1(m,x) + (1.0 - (m/x)^2)*besselj(m,x)*hankelh1(m,x)

d3D(x,m) = x * diffsbesselj(m,x) * (x * diffshankelh1(m,x) + shankelh1(m,x)) + (x^2 + - m * (m+1)) * sbesselj(m,x) * shankelh1(m,x)

function far_field_pattern(ω::T, medium::Acoustic{T,2}, species::Species{2}; basis_order = 2) where T<:Number

    Zs = - get_t_matrices(medium, species, ω, basis_order)
    num_density_inv = one(T)/sum(number_density.(species))

    far_field(θ::T) where T <: Number = -num_density_inv*sum(
        number_density(species[i])*Zs[i][n+basis_order+1,n+basis_order+1]*exp(im*θ*n)
    for i in eachindex(species), n = -basis_order:basis_order)

    return far_field
end

function diff_far_field_pattern(ω::T, medium::Acoustic{T,2}, species::Species{2}; tol=1e-6, basis_order = 2, verbose = false, kws...) where T<:Number

    Zs = - get_t_matrices(medium, species, ω, basis_order)
    num_density_inv = one(T) / sum(number_density.(species))

    far_field(θ::T) where T <: Number = -num_density_inv*sum(
        number_density(species[i])*Zs[i][n+basis_order+1,n+basis_order+1]*im*n*exp(im*θ*n)
    for i in eachindex(species), n=-basis_order:basis_order)

    return far_field
end

function pair_field_pattern(ω::T, medium::Acoustic{T,2}, species::Species{2}; tol::T = T(1e-6), basis_order = 2) where T<:Number

    Zs = - get_t_matrices(medium, species, ω, basis_order)
    num_density_inv = one(T)/sum(number_density.(species))

    pair_field(θ::T) where T <: Number = -T(π)*num_density_inv^(2.0)*sum(
        begin
            a12 = (outer_radius(species[i]) * species[i].separation_ratio + outer_radius(species[j]) * species[j].separation_ratio )
            number_density(species[i]) * number_density(species[j]) * a12^2.0 * d2D(a12*ω/medium.c,m-n) * Zs[i][n+basis_order+1,n+basis_order+1] * Zs[j][m+basis_order+1,m+basis_order+1]*exp(im*m*θ)
        end
    for i in eachindex(species), j in eachindex(species),
    n = -basis_order:basis_order, m = -basis_order:basis_order)

    return pair_field
end
