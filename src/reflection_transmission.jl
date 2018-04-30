"The average reflection coefficient, can return a vector or just one"
function reflection_coefficient(ωs::Union{T,AbstractVector{T}}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number

    k_effs = wavenumber(ωs, medium, species; kws...)
    return reflection_coefficient(ωs, k_effs, medium, species; kws...)
end

"The average reflection coefficients"
reflection_coefficient(ωs::AbstractVector{T}, k_effs::Vector{Complex{T}}, medium::Medium{T}, species::Vector{Specie{T}}; kws...) where T<:Number =
[ reflection_coefficient(ωs[i], k_effs[i], medium, species; kws...) for i in eachindex(ωs)]

"The average reflection coefficient"
function reflection_coefficient(ω::T, k_eff::Complex{T}, medium::Medium{T}, species::Vector{Specie{T}};
        θ_inc::T = 0.0, tol = 1e-8, kws...) where T<:Number

    k = ω/medium.c
    θ_eff = transmission_angle(k, k_eff, θ_inc; tol = tol)
    As = transmission_scattering_coefficients(ω, k_eff, medium, species; tol = tol, θ_inc = θ_inc, kws...)

    θ_ref = pi - θ_eff - θ_inc
    S = length(species)
    ho = Int((size(As,1)-1)/2) # largest hankel order

    R = 2.0im/(k*cos(θ_inc)*(k*cos(θ_inc) + k_eff*cos(θ_eff)))
    R = R*sum(
        exp(im*n*θ_ref)*species[l].num_density*As[n+ho+1,l]*Zn(ω,species[l],medium,n)
    for n=-ho:ho, l=1:S)

    return R
end

function transmission_scattering_coefficients_field(ω::T, k_eff::Complex{T}, medium::Medium{T}, species::Vector{Specie{T}};
        max_hankel_order=3, θ_inc=0.0, kws...) where T<:Number

    As = transmission_scattering_coefficients(ω, k_eff, medium, species;
            θ_inc=θ_inc, hankel_order=max_hankel_order, kws...)
    k = ω/medium.c
    ho = max_hankel_order
    S = length(species)
    θ_eff = transmission_angle(k, k_eff, θ_inc)
    AA(x) = [im^Float64(m)*exp(-im*m*θ_eff)*As[m+ho+1,s]*exp(im*k_eff*cos(θ_eff)*x) for m=-ho:ho, s=1:S]

    return AA
end

"The average transmitted scattering coefficients for a given effective wavenumber k_eff. Assumes there exists only one k_eff.
The function returns an array A, where
AA(x,y,m,s) = im^m*exp(-im*m*θ_eff)*A[m + max_hankel_order +1,s]*exp(im*k_eff*(cos(θ_eff)*x + sin(θ_inc)*y))
where (x,y) are coordinates in the halfspace, m-th hankel order, s-th species,  and AA is the ensemble average scattering coefficient."
function transmission_scattering_coefficients(ω::T, k_eff::Complex{T}, medium::Medium{T}, species::Vector{Specie{T}};
            hankel_order = :auto,
            max_hankel_order = 10,
            radius_multiplier = 1.005,
            tol = 1e-8, θ_inc::T = 0.0,
            kws...) where T<:Number

        k = ω/medium.c
        ho = hankel_order
        S = length(species)
        θ_eff = transmission_angle(k, k_eff, θ_inc; tol = tol)

        as = radius_multiplier*[(s1.r + s2.r) for s1 in species, s2 in species]
        function M(keff,j,l,m,n)
            (n==m ? 1.0:0.0)*(j==l ? 1.0:0.0) + 2.0pi*species[l].num_density*Zn(ω,species[l],medium,n)*
                Nn(n-m,k*as[j,l],keff*as[j,l])/(k^2.0-keff^2.0)
        end

        if hankel_order == :auto
            ho = -1 + sum([ tol .< norm([M(0.9*k + 0.1im,j,l,1,n) for j = 1:S, l = 1:S]) for n=0:max_hankel_order ])
        else ho = hankel_order
        end

        # calculate effective amplitudes
        MM_svd = svd(
            reshape(
                [M(k_eff,j,l,m,n) for m in -ho:ho, j in 1:S, n in -ho:ho, l in 1:S]
            , ((2ho+1)*S, (2ho+1)*S))
        )

        A_null = MM_svd[3][:,(2ho+1)*S] # eignvector of smallest eigenvalue
        A_null = reshape(A_null, (2*ho+1,S)) # A_null[:,j] = [A[-ho,j],A[-ho+1,j],...,A[ho,j]]

        sumAs = 2*sum(
                exp(im*n*(θ_inc - θ_eff))*Zn(ω,species[l],medium,n)*species[l].num_density*A_null[n+ho+1,l]
        for n = -ho:ho, l = 1:S)
        x = im*k*cos(θ_inc)*(k_eff*cos(θ_eff) - k*cos(θ_inc))/sumAs

        return A_null*x
end
