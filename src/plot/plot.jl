@recipe function plot(avg_wave::AverageWave{T};
        hankel_indexes = -avg_wave.hankel_order:avg_wave.hankel_order,
        apply = real) where T<:AbstractFloat

    ho = avg_wave.hankel_order

    for n in hankel_indexes

        apply_field = apply.(avg_wave.amplitudes[:,n+ho+1,1])

        @series begin
            label --> "$apply Hankel = $n"
            (avg_wave.X, apply_field)
        end
    end
end

@recipe function plot(x::AbstractVector{T}, wave_effs::Vector{EffectiveWave{T}};
        hankel_indexes = -wave_effs[1].hankel_order:wave_effs[1].hankel_order,
        apply = real) where T<:AbstractFloat

    ho = wave_effs[1].hankel_order
    avg_wave_effs = [AverageWave(1.0, wave, x) for wave in wave_effs]
    amps = sum(avg_wave_effs[i].amplitudes[:,:,:] for i in eachindex(avg_wave_effs))
    avg_wave_eff = AverageWave(ho, x, amps)

    for n in hankel_indexes

        apply_field = apply.(avg_wave_eff.amplitudes[:,n+ho+1,1])

        @series begin
            label --> "$apply H = $n"
            (x, apply_field)
        end
    end
end
