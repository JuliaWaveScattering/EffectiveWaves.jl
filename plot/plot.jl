include("plot_match.jl")

@recipe function plot(discrete_wave::DiscretePlaneWaveMode{T};
        hankel_indexes = 0:discrete_wave.basis_order,
        apply = real) where T<:AbstractFloat

    ho = discrete_wave.basis_order

    for n in hankel_indexes

        apply_field = apply.(discrete_wave.amplitudes[:,n+ho+1,1])

        @series begin
            label --> "$apply Hankel = $n"
            seriestype --> :scatter
            (discrete_wave.x, apply_field)
        end
    end
end

@recipe function plot(wave_effs::Vector{w}) where w <: AbstractWaveMode
    k_effs = [w.wavenumber for w in wave_effs]
    maxamp = maximum(norm(w.eigenvectors) for w in wave_effs)

    alphas = map(wave_effs) do w
        norm(w.eigenvectors)/maxamp
    end

    @series begin
        ylims --> (0,Inf)
        xguide --> "Re k_eff"
        yguide --> "Im k_eff"
        seriestype --> :scatter
        label --> ""
        markerstrokealpha --> 0.2
        # markercolor --> :blue
        markeralpha --> alphas
        k_effs
    end
end

@recipe function plot(x::AbstractVector{T}, wave_eff::EffectivePlaneWaveMode{T}) where T<:AbstractFloat
    @series begin
        (x, [wave_eff])
    end
end

@recipe function plot(x::AbstractVector, wave_effs::Vector{E};
        halfspace = Halfspace((Dim==2) ? [-one(T),zero(T)] : [zeros(T,Dim-1);-one(T)]),
        hankel_indexes = 0:wave_effs[1].basis_order,
        apply = real) where E<:EffectivePlaneWaveMode{T,Dim} where {T,Dim}

    avg_wave_eff = DiscretePlaneWaveMode(x, wave_effs, halfspace)
    ho = avg_wave_eff.basis_order

    for n in hankel_indexes
        apply_field = apply.(avg_wave_eff.amplitudes[:,n+ho+1,1])

        @series begin
            label --> "$apply Hankel = $n"
            (x, apply_field)
        end
    end
end
