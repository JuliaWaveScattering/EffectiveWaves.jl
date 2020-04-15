# # Make the planar case the default
# eigensystem(ω::T, medium::PhysicalMedium{T,Dim}, species::Species{T,Dim}; kws...) where Dim =  eigensystem(ω, medium, species, PlanarSymmetry{Dim}(); kws...)

eigensystem(ω::T, source::AbstractSource{T}, material::Material; kws...) where T<:AbstractFloat = eigensystem(ω, source.medium, material.species, setupsymmetry(source,material); kws...)

"Calculates the effective wavenumbers and return Vector{EffectivePlaneWaveMode}."
# function effective_wavemodes(ω::T, source::AbstractSource{T}, material::Material{Dim,S,Sps}; kws...) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}} # keeps getting "Unreachable reached" error
function effective_wavemodes(ω::T, source::PlaneSource{T,Dim}, material::Material{Dim,S,Sps}; kws...) where {T,Dim,S<:Shape{T,Dim},Sps<:Species{T,Dim}} # without the parametric types we get a "Unreachable reached" error

    k_effs = wavenumbers(ω, source.medium, material.species; kws... )
    wave_effs = [
        effective_wavemode(ω, k_eff, source, material; kws...)
    for k_eff in k_effs]

    return wave_effs
end
