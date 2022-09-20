struct PercusYevick{Dim} <: PairCorrelationType
    "Relative tolerance used when calculating the Percus-Yevick approximation"
    rtol::Float64
    "Maximum number of quadture evaluations when calculating the Percus-Yevick approximation"
    maxevals::Int
    "Maximum number of points for the pair correlation"
    maxsize::Int
end

struct HoleCorrection <: PairCorrelationType
end

PercusYevick(Dim; rtol::AbstractFloat = 1e-2, maxevals::Int = Int(2e4), maxsize::Int = 50) = PercusYevick{Dim}(rtol, maxevals, maxsize)

"""
    DiscretePairCorrelation

Represents the pair correlation between two types of species, which could be the same.
"""
struct DiscretePairCorrelation <: PairCorrelation
    "distance between particles centres"
    r::Vector{Float64}
    "variation of the pair correlation from 1 (uncorrelated case)"
    dp::Vector{Float64}

    function DiscretePairCorrelation(r::AbstractVector, dp::AbstractVector; tol::AbstractFloat = 1e-3)
        if !isempty(dp) && size(dp) != size(r)
            @error "the size of vector of distances `r` (currently $(size(r))) should be the same as the size of the pair-correlation variation `dp` (currently $(size(dp)))."
        end
        if !isempty(dp) && abs(dp[end]) > tol
            @warn "For the pair-correlation to be accurate, we expect it to be long enough (in terms of the distance `r`) such that the particle positions become uncorrelatd. They become uncorrelated when `dp[end]` tends to zero."
        end
        new(r,dp)
    end
end

PairCorrelation(r::AbstractVector{T},dp::AbstractVector{T}) where T <: AbstractFloat = DiscretePairCorrelation(r,dp)

function DiscretePairCorrelation(s1::Specie, s2::Specie)

    # Have no pair correlation, then only hole correction will be used.

    T = typeof(outer_radius(s1))

    # a12 = outer_radius(s1) * s1.exclusion_distance + outer_radius(s2) * s2.exclusion_distance
    # r = [a12]
    # dp = [zero(typeof(a12))]

    return DiscretePairCorrelation(T[],T[])
end

"""
    DiscretePairCorrelation(s::Specie, pc::PercusYevick, distances::AbstractVector)

Generates a DiscretePairCorrelation for the specie `s` by using the Percus-Yevick approximation. This distribution assumes particles are distributed accoriding to a random uniform distribution, and that particles can not overlap.
"""
function DiscretePairCorrelation(s1::Specie{3}, pc::PairCorrelationType;
        distances::AbstractVector{T} where T<:AbstractFloat = Float64[]
    )

    r1 = outer_radius(s1) * s1.exclusion_distance

    R = 2r1
    numdensity = number_density(s1)

    automatic_dist = if isempty(distances)
        distances = R:(pc.rtol):(10R)
        if length(distances) > pc.maxsize
            distances = distances[1:pc.maxsize]
        end

        true
    else false
    end

    dp = calculate_pair_correlation(R, distances, pc;
        number_density = numdensity
    )

    if automatic_dist
        i = findfirst(reverse(abs.(dp)) .> pc.rtol)
        if isnothing(i)
            dp = typeof(dp)[]
            distances = typeof(dp)[]
        else
            dp = dp[1:end-i+1]
            distances = distances[1:end-i+1]
        end
    end

    return DiscretePairCorrelation(distances, dp; tol = pc.rtol)
end

function DiscretePairCorrelation(s1::Specie{3}, s2::Specie{3}, pc::PairCorrelationType;
        distances::AbstractVector{T} where T<:AbstractFloat = Float64[]
    )

    r1 = outer_radius(s1) * s1.exclusion_distance
    r2 = outer_radius(s2) * s2.exclusion_distance

    if r1 != r2
        @warn "Percus-Yevick approximation has only been implemented for particles of the same size. Will use the average size of both particles and the combine volume fraction"
    end

    R = r1 + r2
    numdensity = number_density(s1) + number_density(s2)

    automatic_dist = if isempty(distances)
        distances = R:(10pc.rtol):(10R)
        true
    else false
    end

    dp = calculate_pair_correlation(R, distances, pc;
        number_density = numdensity
    )

    if automatic_dist
        i = findfirst(reverse(abs.(dp)) .> pc.rtol)
        dp = dp[1:end-i]
        distances = distances[1:end-i]
    end

    return DiscretePairCorrelation(distances, dp; tol = pc.rtol)
end

function calculate_pair_correlation(R::T, distances::AbstractVector{T}, pc::PercusYevick{3};
        number_density::T = 0
    ) where T

    rtol = pc.rtol;
    maxevals = pc.maxevals;

    # Note that f if the volume fraction of the species including the exclusion volume around each particle
    f = number_density * π * R^3 / 6
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

    maxker = ker.(distances) |> maximum

    d = 5 * R
    max_xs = LinRange(d, 30 * d, 200)

    imax = findfirst(ker.(max_xs) ./ maxker .< rtol)
    max_x = isnothing(imax) ? max_xs[end] : max_xs[imax]

    (I,E) = hquadrature(ker, eps(T), max_x;
        rtol = rtol, maxevals = maxevals
    )
    if abs(E/I) > rtol
        @warn "For rtol = $rtol, the estimated relative error when calculating the integral of the Percus-Yevick is: ", abs(E/I)
    end

    dp = map(distances) do r
        if r < R return zero(T) end

        ker = ker_fun(r)
        int_ker = hquadrature(ker, 10eps(T), max_x;
            rtol = rtol, maxevals = maxevals
        )[1]

        η = r / R
        9f * (1 + f) / (2 * η * (1 - f)^3) + 2 / (η*pi) * int_ker
    end

    return dp
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
