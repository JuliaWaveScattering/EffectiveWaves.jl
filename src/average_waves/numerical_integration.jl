
# ints_trap = map(hs) do h
#     x = y1:h:y2
#     -(g(x[1]) + g(x[end]))*h/2.0 + sum(g.(x)*h)
# end

function integration_scheme(x::AbstractVector{T}; scheme = :trapezoidal, kws...) where T<:AbstractFloat
    if length(x) == 1
        return [0.]
    end

    if scheme == :trapezoidal
        trap_scheme(x; kws...)
    elseif scheme == :simpson
        simpson_scheme(x; kws...)
    else
        warn("Integration scheme $scheme unknown, will use trapezoidal")
        trap_scheme(x; kws...)
    end
end

function trap_scheme(x::AbstractVector{T}; x0 = first(x), xn = last(x)) where T<:AbstractFloat

    inds = indices(x,1)

    h = (x[inds[2]]-x[inds[1]])
    σs = x.*zero(T) .+ h
    σs[inds[1]] -= h/T(2)

    σs[end] -= h/T(2)
    # accounts for the ends
    σs[inds[1]] +=   (x[inds[1]] - x0)*(x[inds[2]] - x0 + h)/(2.0h)
    σs[inds[2]] +=  -(x0 - x[inds[1]])^T(2)/(T(2)*h)
    σs[end] +=  (xn - x[end])*(xn - x[end-1] + h)/(2h)
    σs[end-1] += - (xn - x[end])^T(2)/(T(2)*h)

    σs
end

function simpson_scheme(x::AbstractVector{Float64}; x0=x[1], xn=x[end])
    if iseven(length(x))
        warn("Simpson's integration scheme is designed for an odd number of mesh points")
    end

    h = (x[2]-x[1])
    σs = [iseven(j) ? 4.0 : 2.0 for j in eachindex(x)]
    σs = σs*h/3.0

    # accounts for the ends
    σs[1] = h/3.0
    σs[1] +=  (x[1] - x0)*(x[2] - x0 + h)/(2.0h)
    σs[2] +=  - (x0-x[1])^2.0/(2.0h)
    σs[end] = h/3.0
    σs[end] +=  (xn-x[end])*(xn - x[end-1] + h)/(2h)
    σs[end-1] += - (xn-x[end])^2.0/(2h)

    σs
end
