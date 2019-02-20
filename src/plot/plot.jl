include("plot_match.jl")

@recipe function plot(avg_wave::AverageWave{T};
        hankel_indexes = 0:avg_wave.hankel_order,
        apply = real) where T<:AbstractFloat

    ho = avg_wave.hankel_order

    for n in hankel_indexes

        apply_field = apply.(avg_wave.amplitudes[:,n+ho+1,1])

        @series begin
            label --> "$apply Hankel = $n"
            seriestype --> :scatter
            (avg_wave.x, apply_field)
        end
    end
end

@recipe function plot(wave_effs::Vector{EffectiveWave{T}}) where T<:AbstractFloat
    k_effs = [w.k_eff for w in wave_effs]
    maxamp = maximum(norm(w.amplitudes) for w in wave_effs)

    alphas = map(wave_effs) do w
        norm(w.amplitudes)/maxamp
    end

    @series begin
        ylims --> (0,Inf)
        xlab --> "Re k_eff"
        ylab --> "Im k_eff"
        seriestype --> :scatter
        label --> ""
        markerstrokealpha --> 0.2
        # markercolor --> :blue
        markeralpha --> alphas
        k_effs
    end
end

@recipe function plot(x::AbstractVector{T}, wave_eff::EffectiveWave{T}) where T<:AbstractFloat
    @series begin
        (x, [wave_eff])
    end
end

@recipe function plot(x::AbstractVector, wave_effs::Vector{E};
        hankel_indexes = 0:wave_effs[1].hankel_order,
        apply = real) where E<:EffectiveWave

    wave_eff = AverageWave(x, wave_effs)
    ho = wave_eff.hankel_order

    for n in hankel_indexes

        apply_field = apply.(wave_eff.amplitudes[:,n+ho+1,1])

        @series begin
            label --> "$apply Hankel = $n"
            (x, apply_field)
        end
    end
end
