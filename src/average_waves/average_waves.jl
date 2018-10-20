"A type for the ensemble average scattering coefficients.
Here they are discretised in terms of the depth x of the halfspace"
mutable struct AverageWave{T<:AbstractFloat}
    hankel_order::Int # largest hankel order
    x::Vector{T} # spatial mesh
    amplitudes::Array{Complex{T}} # a matrix of the scattering amplitudes, size(A_mat) = (length(x), 2hankel_order +1)
    # Enforce that the dimensions are correct
    function AverageWave{T}(hankel_order::Int, x::Vector{T}, amplitudes::Array{Complex{T}}) where T <: AbstractFloat
        if (length(x), 2*hankel_order+1) != size(amplitudes)[1:2]
            error("The amplitudes of AverageWave does not satisfy size(amplitudes)[1:2] == (length(X), 2*hankel_order+1)")
        end
        new(hankel_order,x,amplitudes)
    end
end

AverageWave(M::Int, x::AbstractVector{T}, as::AbstractArray{Complex{T}}) where T<:AbstractFloat = AverageWave{T}(M,collect(x),collect(as))

function AverageWave(x::AbstractVector{T}, A_mat::Array{Complex{T}}) where T<:Number
    AverageWave(Int((size(A_mat,2)-1)/2), collect(x), A_mat)
end

function average_error(avg_w::AverageWave)
    ddf = circshift(avg_w.amplitudes, (2,0,0)) - 2*circshift(avg_w.amplitudes, (1,0,0)) + avg_w.amplitudes
    h = (avg_w.x[2] - avg_w.x[1])
    ddf = ddf[3:end,:,:]/(h^2)
    max_ddf = maximum(abs.(ddf))

    # trapezoidal error
    return max_ddf*h^3/12
end

"Calculates an AverageWave from one EffectiveWave"
function AverageWave(xs::AbstractVector{T}, wave_eff::EffectiveWave{T}) where T<:Number

    amps = wave_eff.amplitudes
    ho = wave_eff.hankel_order
    θ_eff = wave_eff.θ_eff

    S = size(amps,2)

    average_amps = [
        im^T(m)*exp(-im*m*θ_eff)*amps[m+ho+1,s]*exp(im*wave_eff.k_eff*cos(θ_eff)*x)
    for x in xs, m=-ho:ho, s=1:S]

    return AverageWave(ho,xs,average_amps)
end

"Calculates an AverageWave from a vector of EffectiveWave"
function AverageWave(xs::AbstractVector{T}, wave_effs::Vector{EffectiveWave{T}}) where T<:Number

    ho = wave_effs[1].hankel_order
    avg_wave_effs = [AverageWave(xs, wave) for wave in wave_effs]
    amps = sum(avg_wave_effs[i].amplitudes[:,:,:] for i in eachindex(avg_wave_effs))

    return AverageWave(ho, xs, amps)
end



"Numerically solved the integral equation governing the average wave. Optionally can use wave_eff to approximate the wave away from the boundary."
function AverageWave(ω::T, medium::Medium{T}, specie::Specie{T};
        radius_multiplier::T = 1.005,
        x::AbstractVector{T} = [zero(T)],
        tol::T = T(1e-4),
        wave_effs::Vector{EffectiveWave{T}} = [zero(EffectiveWave{T})],
        max_size::Int = 1000,
        kws...
    ) where T<:Number

    k = real(ω/medium.c)

    if x == [zero(T)]
        if maximum(abs(w.k_eff) for w in wave_effs) == zero(T)
            wave_effs = effective_waves(real(ω/medium.c), medium, [specie];
                radius_multiplier=radius_multiplier, tol=tol, mesh_points=2, kws...)
        end
        # estimate a large coarse non-dimensional mesh based on the lowest attenuating effective wave
        a12 = T(2)*radius_multiplier*specie.r
        x = x_mesh(wave_effs[1]; tol = tol,  a12 = a12, max_size=max_size)
    end

    X = x.*k
    (MM_quad,b_mat) = average_wave_system(ω, X, medium, specie; tol = tol, kws...);

    M = Int( (size(b_mat,2) - 1)/2 )
    J = length(collect(X)) - 1

    len = (J + 1) * (2M + 1)
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    As = MM_mat\b
    As_mat = reshape(As, (J+1, 2M+1, 1))

    return AverageWave(M, collect(X)./k, As_mat)
end

"note that this uses the non-dimensional X = k*depth"
function average_wave_system(ω::T, X::AbstractVector{T}, medium::Medium{T}, specie::Specie{T};
        θin::Float64 = 0.0, tol::T = 1e-6,
        radius_multiplier::T = 1.005,
        scheme::Symbol = :trapezoidal,
        hankel_order::Int = maximum_hankel_order(ω, medium, [specie]; tol = tol),
        t_vecs::Vector{Vector{Complex{T}}} = t_vectors(ω, medium, [specie]; hankel_order = hankel_order),
        kws...
    ) where T<:AbstractFloat

    k = real(ω/medium.c)
    a12k = radius_multiplier*T(2)*real(k*specie.r);
    M = hankel_order;

    J = length(X) - 1
    h = X[2] - X[1]

    PQ_quad = intergrand_kernel(X, a12k; M = M, θin = θin, scheme=scheme);

    MM_quad = [
        (specie.num_density/(k^2))*t_vecs[1][m+M+1]*PQ_quad[l,m+M+1,j,n+M+1] + ( (m==n && j==l) ? 1.0+0.0im : 0.0+0.0im)
    for  l=1:(J+1), m=-M:M, j=1:(J+1), n=-M:M];

    b_mat = [
        -t_vecs[1][m+M+1]*exp(im*X[l]*cos(θin))*exp(im*m*(pi/2.0 - θin))
    for l = 1:(J+1), m = -M:M]

    return (MM_quad,b_mat)
end

"Returns x the mesh used to discretise the integral equations."
function x_mesh(wave_eff_long::EffectiveWave{T}, wave_eff_short::EffectiveWave{T} = wave_eff_long;
        tol::T = T(1e-5),  a12::T = T(Inf),
        max_size::Int = 1000,
        min_size::Int = 5,
        max_x::T = (-log(tol))/abs(cos(wave_eff_long.θ_eff)*imag(wave_eff_long.k_eff))
) where T<:AbstractFloat

    #= The default min_X result in:
        abs(exp(im*min_X*cos(θ_effs[end])*k_effs[end]/k)) < tol
    =#
    # estimate a reasonable derivative based on more rapidly varying wave_eff_short.
    df = abs(wave_eff_short.k_eff * cos(wave_eff_short.θ_eff))

    # Based on Simpson's rule
        # dX  = (tol*90 / (df^4))^(1/5)
    # Based on trapezoidal integration
        dx  = (tol * T(30) / (df^2))^(1/3)

    # prioritise a mesh fine enough to give accurate discrete integrals
    if dx > a12 dx = a12 end

    if max_x/dx + 1 > max_size
        if dx == a12
            max_x = (max_size - 1)*dx
            @warn("The mesh max_size = $max_size which was too small for tol = $tol. Will shrink meshed region.")
        else # otherwise priortise shrink max_x and increase dx proportionally
            a = sqrt(max_x/((max_size - 1)*dx))
            dx = min(dx*a,a12)
            max_x = (max_size - 1)*dx
            @warn("The mesh max_size = $max_size which was too small for tol = $tol. Will make a smaller and coarser mesh.")
        end
    elseif max_x/dx + 1 < min_size
        max_x = (min_size - 1)*dx
    end

    # if whole correction length a12k was given, then make dX/a12k = integer
    if a12 != T(Inf) && dx < a12
        n = ceil(a12 / dx)
        dx = a12/n
    end

    return zero(T):dx:max_x
end
