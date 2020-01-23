"""
    EffectivePlaneWaveMode(medium::P, amplitude::SVector, wavevector::SVector)

Is a struct that represents one mode of the effective scattering coefficients for plane wave symmetry.

This wave mode has frequency ω and has has the value ``A e^{i \\mathbf v \\cdot \\mathbf x }`` at the point ``\\mathbf x``, where ``A`` are the amplitudes and ``\\mathbf v`` is the wavevector.
"""
struct EffectivePlaneWaveMode{T<:AbstractFloat,Dim}
    ω::T
    basis_order::Int
    amplitudes::Array{Complex{T}} # the effective ampliudes
    wavevector::SVector{Dim,Complex{T}} # the effective wavevector
end

function EffectivePlaneWaveMode(ω::T, amps::Array{Complex{T}}, wavevector::Array{Complex{T}}) where T<:AbstractFloat
    EffectivePlaneWaveMode(ω,Int( (size(amps,1) - 1) / 2 ), amps, wavevector)
end

zero(W::Type{EffectivePlaneWaveMode{T}}) where {T<:AbstractFloat} = EffectivePlaneWaveMode(zero(T),0,[zero(Complex{T})],zero(Complex{T}),zero(Complex{T}))

# effective_wavemodes(ω::T, medium::PhysicalMedium{T}, specie::Specie{T}; kws...) where T<:AbstractFloat =  effective_wavemodes(ω, medium, [specie]; kws...)

"Calculates the effective wavenumbers and return Vector{EffectivePlaneWaveMode}."
function effective_wavemodes(ω::T, source::AbstractSource, material::Material; kws...) where {T<:AbstractFloat,Dim}

    k_effs = wavenumbers(ω, source.medium, material.species; kws... )
    wave_effs = [
        EffectivePlaneWaveMode(ω, k_eff, psource, material; kws...)
    for k_eff in k_effs]

    return wave_effs
end
