"A type for the ensemble average scattering coefficients.
Here they are discretised in terms of the depth x of the halfspace"
mutable struct DiscretePlaneWaveMode{T<:AbstractFloat}
    basis_order::Int # largest hankel order
    x::Vector{T} # spatial mesh
    amplitudes::Array{Complex{T}} # a matrix of the scattering amplitudes, size(A_mat) = (length(x), 2basis_order +1, number_of_species)
    # Enforce that the dimensions are correct
    function DiscretePlaneWaveMode{T}(basis_order::Int, x::Vector{T}, amplitudes::Array{Complex{T}}) where T <: AbstractFloat
        if (length(x), 2*basis_order+1) != size(amplitudes)[1:2]
            error("The amplitudes of DiscretePlaneWaveMode does not satisfy size(amplitudes)[1:2] == (length(X), 2*basis_order+1)")
        end
        new(basis_order,x,amplitudes)
    end
end

DiscretePlaneWaveMode(M::Int, x::AbstractVector{T}, as::AbstractArray{Complex{T}}) where T<:AbstractFloat = DiscretePlaneWaveMode{T}(M,collect(x),collect(as))

function DiscretePlaneWaveMode(x::AbstractVector{T}, A_mat::Array{Complex{T}}) where T<:Number
    DiscretePlaneWaveMode(Int((size(A_mat,2)-1)/2), collect(x), A_mat)
end

"Approximates the error in DiscretePlaneWaveMode.amplitudes due to the mesh DiscretePlaneWaveMode.x."
function discretewave_error(avg_w::DiscretePlaneWaveMode)
    ddf = circshift(avg_w.amplitudes, (2,0,0)) - 2*circshift(avg_w.amplitudes, (1,0,0)) + avg_w.amplitudes
    h = (avg_w.x[2] - avg_w.x[1])
    ddf = ddf[3:end,:,:]/(h^2)
    max_ddf = maximum(abs.(ddf))

    # trapezoidal error
    return max_ddf*h^3/12
end

"Calculates an DiscretePlaneWaveMode from one EffectiveWave"
function DiscretePlaneWaveMode(xs::AbstractVector{T}, wave_eff::EffectivePlaneWaveMode{T}, halfspace::Halfspace{T}) where T<:Number

    amps = wave_eff.eigenvectors
    ho = wave_eff.basis_order

    θ_eff = transmission_angle(wave_eff, halfspace)
    kcos_eff = wave_eff.wavenumber * dot(-conj(halfspace.normal), wave_eff.direction)

    S = size(amps,2)
    P = size(amps,3)
    if P > 1
        @warn "The plane-wave has more than one eigenvector per wavenumber. This case has not been implemented so using only the first. This occured for: ω = $(wave_eff.ω) and wavenumber = $(wave_eff.wavenumber)"
    end

    average_amps = [
        im^T(m)*exp(-im*m*θ_eff) * amps[m+ho+1,s,1] * exp(im*kcos_eff*x)
    for x in xs, m=-ho:ho, s=1:S]

    return DiscretePlaneWaveMode(ho,xs,average_amps)
end

"Calculates an DiscretePlaneWaveMode from a vector of EffectiveWave"
function DiscretePlaneWaveMode(xs::AbstractVector{T}, wave_effs::Vector{EffectivePlaneWaveMode{T,Dim}}, halfspace::Halfspace{T}) where {T<:Number,Dim}

    N = wave_effs[1].basis_order
    if !isempty(findall([w.basis_order != N for w in wave_effs]))
        @error "The effective waves $wave_effs have different orders, so they can not be combined into a discrete wave"
    end

    discrete_waves = [DiscretePlaneWaveMode(xs, wave, halfspace) for wave in wave_effs]
    amps = sum(discrete_waves[i].amplitudes[:,:,:] for i in eachindex(discrete_waves))

    return DiscretePlaneWaveMode(N, xs, amps)
end


"Numerically solve the integral equation governing ensemble average waves. Optionally can use wave_eff to approximate the wave away from the boundary."
function DiscretePlaneWaveMode(ω::T, source::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{2,Halfspace{T,2}};
        x::AbstractVector{T} = [zero(T)],
        tol::T = T(1e-4),
        wave_effs::Vector{EffectivePlaneWaveMode{T,2}} = EffectivePlaneWaveMode{T,2}[],
        max_size::Int = 1000,
        kws...
    ) where T<:Number

    k = real(ω/source.medium.c)
    specie = material.microstructure.species[1]

    if x == [zero(T)]
        if isempty(wave_effs)
            wave_effs = WaveModes(ω, source, material; tol=tol, mesh_points=2, kws...)
        end
        # estimate a large coarse non-dimensional mesh based on the lowest attenuating effective wave
        a12 = T(2) * specie.seperation_ratio * outer_radius(specie)
        x = x_mesh(wave_effs[1]; tol = tol,  a12 = a12, max_size=max_size)
    end

    X = x.*k
    (MM_quad,b_mat) = discrete_wave_system(ω, X, source, material; tol = tol, kws...);

    M = Int( (size(b_mat,2) - 1)/2 )
    J = length(collect(X)) - 1

    len = (J + 1) * (2M + 1)
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    As = MM_mat\b
    As_mat = reshape(As, (J+1, 2M+1, 1))

    return DiscretePlaneWaveMode(M, collect(X)./k, As_mat)
end
