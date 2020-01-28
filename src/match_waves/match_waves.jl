"A type for matched waves."
mutable struct MatchPlaneWaveMode{T<:AbstractFloat,Dim}
    effective_wavemodes::Vector{EffectivePlaneWaveMode{T,Dim}}
    discrete_wave::DiscretePlaneWaveMode{T}
    x_match::Vector{T} # waves are matched between discrete_wave.x_match
end

"Calculates the difference between the match of MatchPlaneWaveMode.effective_wavemodes and MatchPlaneWaveMode.discrete_wave. This can be used as a proxi for convergence. "
function match_error(m_wave::MatchPlaneWaveMode{T}, shape::Shape{T}; apply_norm::Function=norm) where T<:AbstractFloat
    avg_eff = DiscretePlaneWaveMode(m_wave.x_match, m_wave.effective_wavemodes,shape)
    j0 = findmin(abs.(m_wave.discrete_wave.x .- m_wave.x_match[1]))[2]
    len = length(m_wave.x_match)

    return apply_norm(m_wave.discrete_wave.amplitudes[j0:end,:,:][:] - avg_eff.amplitudes[:])/len
end

function MatchPlaneWaveMode(ω::T, source::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{2,Halfspace{T,2}};
        tol::T = T(1e-5),
        max_size::Int = 200,
        wave_effs::Vector{EffectivePlaneWaveMode{T,2}} = EffectivePlaneWaveMode{T,2}[],
        x::AbstractVector{T} = [-one(T)],
        L_match::Int = 0,
        kws...
    ) where T<:Number

    θin = transmission_angle(source,material)
    k = real(ω / source.medium.c)

    if isempty(wave_effs)
        wave_effs = effective_wavemodes(k, source, material;
            extinction_rescale = false,
            tol = T(10)*tol,
            kws...)
    end
    basis_order = wave_effs[1].basis_order

    # use non-dimensional effective waves
    wave_non_effs = map(wave_effs) do w
        EffectivePlaneWaveMode(ω,w.basis_order,w.wavevector ./ k, w.amplitudes)
    end

    if first(x) == - one(T)
        # using non-dimensional wave_non_effs and a12k results in non-dimensional mesh X
        a12k = T(2) * k * material.species[1].exclusion_distance * outer_radius(material.species[1])
        L_match, X =  x_mesh_match(wave_non_effs; a12 = a12k, tol = tol, max_size=max_size)
    else
        X = x .* real(k)
        if L_match == 0
            L_match = Int(round(length(X)/2))
        end
    end

    avg_wave_effs = [
        DiscretePlaneWaveMode(X[L_match:L_match+1], w, material.shape)
    for w in wave_non_effs]

    wave_non_effs = map(eachindex(wave_non_effs)) do i
        w = wave_non_effs[i]
        amps = w.amplitudes / norm(avg_wave_effs[i].amplitudes[1,:,1])
        EffectivePlaneWaveMode(ω,w.basis_order,w.wavevector, amps)
    end

    J = length(collect(X)) - 1
    len = (J + 1)  * (2basis_order + 1)

    (MM_quad,b_mat) = discrete_wave_system(ω, X, source.medium, material.species[1]; tol=tol, basis_order=basis_order, θin=θin,  kws...);
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    (LT_mat, ER_mat, b_eff) = match_arrays(ω, wave_non_effs, L_match, X, source, material);

    B = b - ER_mat*b_eff
    As = (ER_mat*LT_mat + MM_mat)\B
    As_mat = reshape(As, (J+1, 2basis_order+1, 1))

    αs = LT_mat*As + b_eff
    # use these αs to correct the magnitude of the amplitudes of the effective waves
    # and re-dimensionalise the effective wavenumbers
    wave_non_effs = map(eachindex(wave_effs)) do i
        w = wave_non_effs[i]
        amps = αs[i] .* w.amplitudes
        vec = k .* w.wavevector
        EffectivePlaneWaveMode(ω,w.basis_order,vec, amps)
    end

    # return MatchPlaneWaveMode(wave_effs, DiscretePlaneWaveMode(basis_order, collect(X)./k, As_mat), collect(X[L_match:end])./k)
    return MatchPlaneWaveMode(wave_non_effs, DiscretePlaneWaveMode(basis_order, real.( collect(X) ./ k), As_mat), real.(collect(X[L_match:end]) ./ k))
end

"Returns (x,L), where x[L:end] is the mesh used to match with wave_effs."
function x_mesh_match(wave_effs::Vector{EffectivePlaneWaveMode{T,Dim}}; kws... ) where {T<:AbstractFloat,Dim}
    # wave_effs[end] establishes how long X should be, while wave_effs[1] establishes how fine the mesh should be.

   # If there is only one wave, then it doesn't make sense to extend the mesh until it decays.
   # Instead we choose, arbitrarily, a quarter of the wavelength.
   # Note, having only one wave is very unusual, but tends to happen in the very low frequency limit.
    x = (length(wave_effs) > 1) ?
        x_mesh(wave_effs[end], wave_effs[1]; kws...) :
        x_mesh(wave_effs[1]; max_x = (pi/2) / abs(sqrt(sum(wave_effs[1].wavevector .^2))), kws...)
    #NOTE previously used abs(cos(wave_effs[1].θ_eff)*abs(wave_effs[1].k_eff)) instead of abs(sqrt(sum(wave_effs[1].wavevector .^2)))

    x_match = x[end]
    x_max = (length(x) < length(wave_effs)*T(1.5)) ?
         x[end] + T(1.5)*(x[2] - x[1])*length(wave_effs) :
         T(2) * x[end]

    x = 0.0:(x[2] - x[1]):x_max

    L = findmin(abs.(x .- x_match))[2]

    return L, x
end
