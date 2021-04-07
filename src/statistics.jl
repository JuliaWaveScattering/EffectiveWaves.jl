function hole_correction_pair_correlation(x1::AbstractVector{T},s1::Specie{T}, x2::AbstractVector{T},s2::Specie{T}) where T
    overlapping = norm(x1 - x2) > outer_radius(s1) * s1.exclusion_distance + outer_radius(s2) * s2.exclusion_distance
    
    return  overlapping ? one(T) : zero(T)
end
