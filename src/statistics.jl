struct PercusYevick <: PairCorrelationType
    # volume_fraction::AbstractFloat
end

struct HoleCorrection <: PairCorrelationType
    # volume_fraction::AbstractFloat
end

"""
    DiscretePairCorrelation

Represents the pair correlation between two types of species, which could be the same.
"""
struct DiscretePairCorrelation <: PairCorrelation
    "distance between particles centres"
    r::AbstractVector{T} where T <: Number
    "variation of the pair correlation from 1 (uncorrelated case)"
    dp::AbstractVector{T} where T <: Number

    function DiscretePairCorrelation(r::AbstractVector,dp::AbstractVector; tol::AbstractFloat = 1e-4)
        if size(dp) != size(r)
            @error "the size of vector of distances `r` (currently $(size(r))) should be the same as the size of the pair-correlation variation `dp` (currently $(size(dp)))."
        end
        if abs(dp[end]) > tol
            @warn "For the pair-correlation to be accurate, we expect it to be long enough (in terms of the distance `r`) such that the particle positions become uncorrelatd. They become uncorrelated when `dp[end]` tends to zero."
        end
        new(r,dp)
    end
end

PairCorrelation(r::AbstractVector{T},dp::AbstractVector{T}) where T <: AbstractFloat = DiscretePairCorrelation(r,dp)

function DiscretePairCorrelation(s::Specie{Dim}, pc::PercusYevick, distances::AbstractVector{T};
        rtol::T = 1e-3, maxevals::Int = Int(2e4)
    ) where {T<:AbstractFloat, Dim}

    R = 2 * outer_radius(s) * s.exclusion_distance

    # Note that f if the volume fraction of the species including the exclusion volume around each particle
    f = number_density(s) * π * R^3 / 6
    α = - (1 + 2f)^2 / (1 - f)^4
    β = 6f * (1 + f/2)^2 / (1 - f)^4
    δ = - f * (1 + 2f)^2 / (2 * (1 - f)^4)

    function F(x)
        A = 24δ / x^6 - 2β / x^4
        B = (α + 2β + 4δ) / x^3 - 24δ / x^5
        C = - (α + β + δ) / x^2 + (2β + 12δ) / x^4 - 24δ / x^6

        return A + B * sin(x) + C * cos(x)
    end

    G(x) = (α + 2β + 4δ) * sin(x) / x^2 - (α + β + δ) * cos(x) / x

    ker_fun(r) = x -> (x * F(x) / (1 - 24f * F(x)) - G(x)) * sin(r * x / R)

    ## Below we estimate the domain of integration needed

    ker = ker_fun(R)

    data = ker.(distances)
    maxker = maximum(data)

    d = 5 * R
    max_xs = LinRange(d, 30 * d, 200)

    imax = findfirst(ker.(max_xs) ./ maxker .< rtol)
    max_x = max_xs[imax]

    (I,E) = hquadrature(ker, eps(T), max_x;
        rtol=rtol, maxevals=maxevals
    )
    println("I = ", maximum(abs.(I)) )
    println("For rtol = $rtol, the estimated relative error when calculating the integral of the Percus-Yevick is: ", abs(E/I))

    dp = map(distances) do r
        ker = ker_fun(r)
        int_ker = hquadrature(ker, 10eps(T), max_x;
            rtol = rtol, maxevals = maxevals
        )[1]

        η = r / R
        9f * (1 + f) / (2 * η * (1 - f)^3) + 2 / (η*pi) * int_ker
    end

    return DiscretePairCorrelation(distances,dp; tol = rtol)
end


function hole_correction_pair_correlation(x1::AbstractVector{T},s1::Specie, x2::AbstractVector{T},s2::Specie) where T <: Number
    overlapping = norm(x1 - x2) > outer_radius(s1) * s1.exclusion_distance + outer_radius(s2) * s2.exclusion_distance

    return  overlapping ? one(T) : zero(T)
end

function smooth_pair_corr_distance(pair_corr_distance::Function, a12::T; smoothing::T = T(0), max_distance::T = T(20*a12),
        polynomial_order::Int = 15, mesh_size::Int = 10*polynomial_order + 4
    ) where T

    if smoothing > T(1.0)
        @warn "smoothing should be in [0,1], setting smoothing = 1"
        smoothing = T(1.0)
    elseif smoothing < T(0.0)
        @warn "smoothing should be in [0,1], setting smoothing = 0"
        smoothing < T(0.0)
    end

    zs = collect(LinRange(zero(T), max_distance, mesh_size))

    inds = findall(abs.(zs .- a12) .< a12*smoothing/T(2) )
    deleteat!(zs,inds)
    data = pair_corr_distance.(zs)

    P = Legendre()
    ls = 0:polynomial_order

    Pmat = P[2zs ./ max_distance .- T(1.0), ls .+ 1];

    # projector_mat = inv(transpose(Pmat) * (Pmat)) * transpose(Pmat);
    # pls = projector_mat * data
    # Pmat * pls ~ data
    pls = Pmat \ data

    return function (z)
        Ps = P[2z / max_distance - T(1.0), ls .+ 1]
        return sum(Ps .* pls)
    end
end

# function constant_number_density_fun(material::Material)
#
#     @warn "the constant_number_density doesn't currently verify if the whole particle is inside the shape"
#
#     prob = function (x1::AbstractVector, s1::Specie)
#         if x1 ∈ sh
#             number_density(s1)
#         else
#             zero(typeof(V))
#         end
#     end
#
#     return prob
# end
"""
    gls_pair_radial_fun(pair_corr_distance::Function, a12::T; polynomial_order::Int, mesh_size::Int)

Return a function ``gls_fun``. For any radial distances ``r_1`` and ``r_2`` we have ``gls = gls_fun(r1,r2)``  such that ``g(r_1,r_2,\\cos \\theta_{12}) = \\sum_{\\ell_1 =0} \\frac{2\\ell_1 + 1}{4\\pi} gls[\\ells+1] P_{\\ell_1}(\\cos \\theta_{12}) ``, where ``g(r_1,r_2,\\cos \\theta_{12})`` is the radially symmetric pair-correlation, so it depends only on the radial distances ``r_1`` and ``r_2``, and the angle between two position vectors ``\\theta_{12}``.

The function `gls_fun` is calculated from the function `pair_corr_distance`, where `pair_corr_distance(sqrt(r1^2 + r2^2 - 2r1 * r2 * cos(θ12)))` gives the pair correlation.
"""
function gls_pair_radial_fun(pair_corr_distance::Function, a12::T;
            polynomial_order::Int = 15,
            mesh_size::Int = 10*polynomial_order + 4,
            sigma_approximation = true
        ) where T

    P = Legendre()
    ls = 0:polynomial_order

    if sigma_approximation
        sigmas = [one(T); sin.(pi .* ls[2:end] ./ (polynomial_order+1)) ./ (pi .* ls[2:end] ./ (polynomial_order+1))]
    else
        sigmas = ones(T,polynomial_order+1)
    end
    S = diagm(sigmas)

    us = LinRange(-1.0, 1.0, mesh_size)
    Pmat = P[us, ls .+ 1];
    Pmat = [Pmat[i] * T(2i[2] - 1) / (4pi) for i in CartesianIndices(Pmat)];

    # Pmat = S * Pmat;

    # projector_mat =  S * inv(transpose(Pmat) * (Pmat)) * transpose(Pmat);

    return function (r1,r2)
        if (r1 + r2 < a12)
            return zeros(typeof(a12), polynomial_order + 1)
        else
            data = pair_corr_distance.(sqrt.(r1^2 .+ r2^2 .- 2r1 .* r2 .* us))
            # pls = projector_mat * data
            # data ~ Pmat * pls
            pls = Pmat \ data

            return S * pls
        end
    end
end

"""
    pair_radial_fun(pair_corr::Function, a12::T; polynomial_order::Int, mesh_size::Int)

Return a function `pair_radial` such that `pair_radial(r1,r2, cos(θ12))` gives the pair correlation particles at the radial distances r1 and r2, with and angle of `θ12` between them.

The function `pair_radial` is calculated from the function `pair_corr_distance`, where `pair_corr_distance(sqrt(r1^2 + r2^2 - 2r1 * r2 * cos(θ12)))` gives the pair correlation.
"""
function pair_radial_fun(pair_corr_distance::Function, a12::T; polynomial_order::Int = 15, kws...) where T

    gls_fun = gls_pair_radial_fun(pair_corr_distance, a12; polynomial_order = polynomial_order, kws...)
    P = Legendre()

    return function (r1,r2,u)
        Pus = P[u, 1:(polynomial_order + 1)] .* (2 .* (0:polynomial_order) .+ 1) ./ (4pi)

        return sum(Pus .* gls_fun(r1,r2))
    end
end

function pair_radial_to_pair_corr(pair_radial::Function)
    function (x1,s1,x2,s2)
        if norm(x1) < 1e-12 || norm(x2) < 1e-12
            pair_radial(norm(x1),norm(x2),0.0)
        else
            pair_radial(norm(x1),norm(x2),dot(x1,x2) / (norm(x1)*norm(x2)))
        end
    end
end
