@recipe function plot(match_wave::MatchPlaneWaveMode{T}) where T <: AbstractFloat

    dx = match_wave.x_match[2] - match_wave.x_match[1]
    max_x = match_wave.x_match[end]
    x = [match_wave.x_match; (max_x:dx:(2*max_x)) .+ dx]
    @series (x, match_wave)
end

@recipe function plot(x::AbstractVector{T}, match_wave::MatchPlaneWaveMode{T};
        basis_order = match_wave.PlaneWaveModes[1].basis_order,
        hankel_indexes = 0:basis_order,
        apply = real, match_region = true) where T <: AbstractFloat

    ho = basis_order
    i1 = findmin(abs.(match_wave.discrete_wave.x .- match_wave.x_match[1]))[2]
    i2 = findmin(abs.(match_wave.discrete_wave.x .- match_wave.x_match[end]))[2]

    max_amp = maximum(apply.(match_wave.discrete_wave.amplitudes[i1:i2,(hankel_indexes) .+ (ho+1),:]))
    min_amp = minimum(apply.(match_wave.discrete_wave.amplitudes[i1:i2,(hankel_indexes) .+ (ho+1),:]))

    max_amp = (max_amp > 0) ? 1.1*max_amp : 0.0
    min_amp = (min_amp < 0) ? 1.1*min_amp : 0.0

    if match_region == true
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
    end

    @series begin
        apply --> apply
        hankel_indexes --> hankel_indexes
        # markeralpha --> 0.2
        markersize --> 4.0
        match_wave.discrete_wave
    end

    @series begin
        apply --> apply
        hankel_indexes --> hankel_indexes
        (x, match_wave.PlaneWaveModes)
    end

end
