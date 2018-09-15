@recipe function plot(avg_wave::AverageWave{T};
        hankel_indexes = -avg_wave.hankel_order:avg_wave.hankel_order,
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
        # ylims --> (0,Inf)
        xlab --> "Re k_eff"
        ylab --> "Im k_eff"
        seriestype --> :scatter
        label --> "wavenumber k_eff"
        markercolor --> :blue
        markerstrokealpha --> alphas
        (real.(k_effs), imag.(k_effs))
    end
end

@recipe function plot(x::AbstractVector{T}, wave_eff::EffectiveWave{T}) where T<:AbstractFloat
    @series begin
        (x, [wave_eff])
    end
end

@recipe function plot(x::AbstractVector{T}, wave_effs::Vector{EffectiveWave{T}};
        hankel_indexes = -wave_effs[1].hankel_order:wave_effs[1].hankel_order,
        apply = real) where T<:AbstractFloat

    wave_eff = AverageWave(x, wave_effs)
    ho = wave_eff.hankel_order

    for n in hankel_indexes

        apply_field = apply.(wave_eff.amplitudes[:,n+ho+1,1])

        @series begin
            label --> "$apply H = $n"
            (x, apply_field)
        end
    end
end

@recipe function plot(match_wave::MatchWave{T}) where T <: AbstractFloat

    dx = match_wave.x_match[2] - match_wave.x_match[1]
    max_x = match_wave.x_match[end]
    x = [match_wave.x_match; (max_x:dx:(2*max_x)) + dx]
    @series (x, match_wave)
end

@recipe function plot(x::AbstractVector{T}, match_wave::MatchWave{T}; hankel_order = match_wave.effective_waves[1].hankel_order,
        hankel_indexes = -hankel_order:hankel_order,
        apply = real) where T <: AbstractFloat

    ho = maximum(hankel_indexes)
    wave_eff = AverageWave(match_wave.x_match, match_wave.effective_waves)
    max_amp = maximum(apply.(wave_eff.amplitudes[:,(hankel_indexes) .+ (ho+1),:]))
    min_amp = minimum(apply.(wave_eff.amplitudes[:,(hankel_indexes) .+ (ho+1),:]))

    max_amp = (max_amp > 0) ? 1.1*max_amp : 0.0
    min_amp = (min_amp < 0) ? 1.1*min_amp : 0.0

    @series begin
        label --> "match region"
        fillalpha --> 0.3
        fill --> (0,:orange)
        line --> 0
        (x1 -> x1, y->max_amp, match_wave.x_match[1],  match_wave.x_match[end])
    end
    @series begin
        label --> ""
        fillalpha --> 0.3
        fill --> (0,:orange)
        line --> 0
        (x1 -> x1, y->min_amp, match_wave.x_match[1],  match_wave.x_match[end])
    end

    @series begin
        apply --> apply
        hankel_indexes --> hankel_indexes
        match_wave.average_wave
    end

    @series begin
        apply --> apply
        hankel_indexes --> hankel_indexes
        (x, match_wave.effective_waves)
    end

end
