reflection_coefficient(ω::T, wavenumber::Complex{T}, psource::PlaneSource{T,Dim}, material::Material{Halfspace{T,Dim}}; kws...) where {T,Dim} = reflection_coefficient(ω, WaveMode(ω, wavenumber, psource, material), psource, material; kws...)

"""
    reflection_coefficient(ω::T, wave_eff::EffectivePlaneWaveMode, psource::PlaneSource{T,2,1,Acoustic}, material::Material{Halfspace}; [x::T = zero(T)])

The reflection coefficient in 2D for acoustics for just one [`EffectivePlaneWaveMode`](@ref)"
"""
function reflection_coefficient(ω::T, wave_eff::EffectivePlaneWaveMode{T}, psource::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{Halfspace{T,2}};
        x::T = zero(T), kws...) where T<:AbstractFloat

    if psource.medium != material.microstructure.medium @error mismatched_medium end

    k = ω / psource.medium.c
    θ_eff = transmission_angle(wave_eff,material)
    θin = transmission_angle(psource,material)

    θ_ref = pi - θ_eff - θin
    S = length(material.microstructure.species)
    ho = wave_eff.basis_order

    kcos_eff = wave_eff.wavenumber * dot(- conj(material.shape.normal), wave_eff.direction)
    kcos_in = k * dot(- conj(material.shape.normal), psource.direction)

    kθ = kcos_in + kcos_eff
    R = 2.0im / (kcos_in * kθ)
    R = R*sum(
        exp(im*n*θ_ref + im*x*kθ) * number_density(material.microstructure.species[l]) *
        # sum(wave_eff.eigenvectors[n+ho+1,l,:])
        wave_eff.eigenvectors[n+ho+1,l,1]
    for n=-ho:ho, l=1:S)

    return R
end

function reflection_coefficient(wavemode::EffectivePlaneWaveMode{T,Dim}, psource::PlaneSource{T,3,1,Acoustic{T,3}}, material::Material{Halfspace{T,3}}) where {T<:AbstractFloat,Dim}


    if psource.medium == material.microstructure.medium

        # Unpacking parameters
        k = wavemode.ω / psource.medium.c

        species = material.microstructure.species
        S = length(species)
        rs = outer_radius.(species)

        basis_order = wavemode.basis_order
        ls, ms = spherical_harmonics_indices(basis_order)

        # make the normal outward pointing
        plate = material.shape
        n = plate.normal / norm(plate.normal)
        if real(dot(n,psource.direction)) > 0
           n = -n
        end

        rθφ = cartesian_to_radial_coordinates(psource.direction)
        Yrefs = spherical_harmonics(basis_order, rθφ[2], rθφ[3]);

        Z0 = dot(-n,plate.origin)
        kcos_in = k * dot(- conj(n), psource.direction)

        kcos_eff = wavemode.wavenumber * dot(- conj(n), wavemode.direction)

        Ramp = - sum(number_density(species[i[2]]) * wavemode.eigenvectors[i] * 2pi * (-1.0)^ms[i[1]] * (1.0im)^(ls[i[1]]-1) *
            Yrefs[i[1]] * exp(im*(kcos_eff + kcos_in)*(Z0 + rs[i[2]])) / ((kcos_in + kcos_eff) * k * kcos_in)
        for i in CartesianIndices(wavemode.eigenvectors))
        Ramp = Ramp * exp(im * kcos_in * Z0)

        return Ramp
    
    else

        # Unpacking parameters
        ρ = psource.medium.ρ
        ρ0 = material.microstructure.medium.ρ
        c = psource.medium.c
        c0 = material.microstructure.medium.c
        ω = wavemode.ω
        k = ω / c
        k0 = ω / c0

        species = material.microstructure.species
        S = length(species)
        rs = outer_radius.(species)

        basis_order = wavemode.basis_order
        ls, ms = spherical_harmonics_indices(basis_order)

        # make the normal outward pointing
        halfspace = material.shape
        n = halfspace.normal / norm(halfspace.normal)
        if real(dot(n,psource.direction)) > 0
            n = -n
        end

        Z0 = dot(-n, halfspace.origin)

        kz = k * dot(-conj(n), psource.direction)
        k0z = sqrt(k0^2 - (k^2 - kz^2))
        k_effz = wavemode.wavenumber * dot(-conj(n), wavemode.direction)

        ζR = (ρ0 * kz - ρ * k0z) / (ρ0 * kz + ρ * k0z)
        ζT = 2ρ * k0z / (ρ0 * kz + ρ * k0z)

        direction_k0 = real(k) * (psource.direction - dot(-conj(n), psource.direction) * psource.direction)
        direction_k0 = direction_k0 + real(k0z) * dot(-conj(n), psource.direction) * psource.direction
        direction_k0 = direction_k0 / norm(direction_k0)

        rθφ = cartesian_to_radial_coordinates(direction_k0)
        Yrefs = spherical_harmonics(basis_order, rθφ[2], rθφ[3]);

        B = - sum(number_density(species[i[2]]) * wavemode.eigenvectors[i] * 2pi * (1.0im)^(ls[i[1]]-1) *
             Yrefs[i[1]] * exp(im * (k_effz + k0z) * (Z0 + rs[i[2]])) / (k0 * k0z * (k_effz + k0z))
        for i in CartesianIndices(wavemode.eigenvectors))

        Ramp = (ζR * field(psource, zeros(T, 3), ω) + ζT * B) * exp(im * kz * Z0)

        return Ramp
    end
end

function reflection_transmission_coefficients(wavemodes::Vector{E}, psource::PlaneSource{T,3,1,Acoustic{T,3}}, material::Material{Plate{T,3}}) where {T<:AbstractFloat,Dim, E<:EffectivePlaneWaveMode{T,Dim}}

    if psource.medium != material.microstructure.medium

        # Unpacking parameters
        ρ = psource.medium.ρ
        ρ0 = material.microstructure.medium.ρ
        c = psource.medium.c
        c0 = material.microstructure.medium.c
        ω = wavemodes[1].ω
        k = ω / c
        k0 = ω / c0
        k_eff = wavemodes[1].wavenumber
        Z = material.shape.width
        G = field(psource, zeros(T, 3), ω)

        species = material.microstructure.species
        S = length(species)
        rs = outer_radius.(species)
        nf = number_density.(species)

        basis_order = maximum(w.basis_order for w in wavemodes)
        if basis_order != minimum(w.basis_order for w in wavemodes)
            @error "expected wavemodes to have the same basis_order"
        end

        ls, ms = spherical_harmonics_indices(basis_order)

        # make the normal outward pointing
        plate = material.shape
        n = plate.normal / norm(plate.normal)
        if real(dot(n,psource.direction)) > 0
            n = -n
        end

        rθφ = cartesian_to_radial_coordinates(psource.direction)
        Ys = spherical_harmonics(basis_order, rθφ[2], rθφ[3])
        lm_to_n = lm_to_spherical_harmonic_index

        Z0 = dot(-n,plate.origin)
        Z1 = Z0 - Z / 2
        Z2 = Z0 + Z / 2

        kz = k * dot(- conj(n), psource.direction)
        k0z = sqrt(k0^2 - k^2 + kz^2)

        # Needed coefficients
        Δ = 2 * k * k0 * ρ * ρ0 * cos(k0 * Z) - 1im * ((k0 * ρ)^2 + (k * ρ0)^2 * sin(k0 * Z))
        D_p = k0 * ρ * (k0 * ρ + k * ρ0) / Δ
        D_m = k0 * ρ * (k0 * ρ - k * ρ0) / Δ
        D_0 = 2 * k * k0 * ρ * ρ0 / Δ
        D_1 = ((k0 * ρ)^2 - (k * ρ0)^2) / 2Δ

        Pr(x::Complex{T}, y::Complex{T}, r::T) = exp(2im * (x + y) / Z) * sin((x + y) * (Z / 2 - r)) / (x + y)

        Bp(Fp::Array{Complex{T}}, Fm::Array{Complex{T}}) = (4π / (k0 * k0z)) * sum(
            1im^(-T(l)) * Ys[lm_to_n(l, m)] * (Pr(k_eff, -k0, rs[j]) * Fp[lm_to_n(l, m), j, p] +
                                            Pr(-k_eff, -k0, rs[j]) * Fm[lm_to_n(l, m), j, p]) * nf[j]
        for p = 1:size(Fp, 3), l = 0:basis_order for m = -l:l, j in eachindex(species))

        Bm(Fp::Array{Complex{T}}, Fm::Array{Complex{T}}) = (4π / (k0 * k0z)) * sum(
            1im^(T(l)) * Ys[lm_to_n(l, m)] * (Pr(k_eff, k0, rs[j]) * Fp[lm_to_n(l, m), j, p] +
                                         Pr(-k_eff, k0, rs[j]) * Fm[lm_to_n(l, m), j, p]) * nf[j]
        for p = 1:size(Fp, 3), l = 0:basis_order for m = -l:l, j in eachindex(species))

        Ramp = D_p * Bm(wavemodes[1].eigenvectors, wavemodes[2].eigenvectors) * exp(-1im * k0 * Z) +
               D_m * Bp(wavemodes[1].eigenvectors, wavemodes[2].eigenvectors) * exp(1im * k0 * Z) + 
               2im * D_1 * G * sin(k0 * Z)

        Tamp = (D_m * Bm(wavemodes[1].eigenvectors, wavemodes[2].eigenvectors) +
                D_p * Bp(wavemodes[1].eigenvectors, wavemodes[2].eigenvectors) +
                D_0 * G) * exp(-1im * k * Z)

        # direction_ref = psource.direction - 2 * dot(n,psource.direction) * n
        # reflected_wave = PlaneSource(psource.medium; direction = direction_ref, amplitude = Ramp)
        # transmitted_wave = PlaneSource(psource.medium; direction = psource.direction, amplitude = Tamp)

        return [Ramp, Tamp]
    else
        # Unpacking parameters
        ω = wavemodes[1].ω
        k = ω / psource.medium.c
        k0 = ω / psource.medium.c

        species = material.microstructure.species
        S = length(species)
        rs = outer_radius.(species)

        basis_order = maximum(w.basis_order for w in wavemodes)
        if basis_order != minimum(w.basis_order for w in wavemodes)
            @error "expected wavemodes to have the same basis_order"
        end

        ls, ms = spherical_harmonics_indices(basis_order)

        # make the normal outward pointing
        plate = material.shape
        n = plate.normal / norm(plate.normal)
        if real(dot(n, psource.direction)) > 0
            n = -n
        end

        rθφ = cartesian_to_radial_coordinates(psource.direction)
        Ys = spherical_harmonics(basis_order, rθφ[2], rθφ[3])

        Z0 = dot(-n, plate.origin)
        Z1 = Z0 - plate.width / 2
        Z2 = Z0 + plate.width / 2

        kz = k * dot(-conj(n), psource.direction)

        RTs = map(wavemodes) do w
            kpz = w.wavenumber * dot(-conj(n), w.direction)

            Rp = sum(
                number_density(species[i[2]]) * w.eigenvectors[i] * 2pi * (1.0im)^(ls[i[1]] - 1) * (-1.0)^ms[i[1]] *
                Ys[i[1]] * (exp(im * (kpz + kz) * (Z2 - rs[i[2]])) - exp(im * (kpz + kz) * (Z1 + rs[i[2]]))) /
                ((kz + kpz) * k * kz)
                for i in CartesianIndices(w.eigenvectors))

            Tp = sum(
                number_density(species[i[2]]) * w.eigenvectors[i] * 2pi * (-1.0)^ls[i[1]] * (1.0im)^(ls[i[1]] + 1) *
                Ys[i[1]] * (exp(im * (kpz - kz) * (Z2 - rs[i[2]])) - exp(im * (kpz - kz) * (Z1 + rs[i[2]]))) /
                ((kz - kpz) * k * kz)
                for i in CartesianIndices(w.eigenvectors))

            [Rp, Tp]
        end

        Ramp, Tamp = (sum(RTs) + [0.0, 1.0]) .* field(psource, zeros(T, 3), ω)

        Ramp = Ramp * exp(im * kz * Z1)
        Tamp = Tamp * exp(im * kz * Z2)

        # direction_ref = psource.direction - 2 * dot(n,psource.direction) * n
        # reflected_wave = PlaneSource(psource.medium; direction = direction_ref, amplitude = Ramp)
        # transmitted_wave = PlaneSource(psource.medium; direction = psource.direction, amplitude = Tamp)

        return [Ramp, Tamp]
    end
end

# function F0(S,j,l,m,n)
#     (S^T(2) - (k*as[j,l])^T(2)) * (n == m ? T(1) : T(0)) * (j == l ? T(1) : T(0)) +
#     T(2) * as[j,l]^T(2) * pi* number_density(species[l]) *t_vecs[l][m+ho+1] * kernelN2D(n-m,k*as[j,l],S)
# end

# kernelN2D(0,k*a12,Z) = k*a12*diffhankelh1(0,k*a12)*besselj(0,Z) - Z*hankelh1(0,k*a12)*diffbesselj(0,Z)
# Q0(S,j,l,m,n) = F0(S,j,l,m,n) / (S^T(2) - (k*as[j,l])^T(2))

# function F0p(S, maxZ::T = maxZ, num_coefs::Int = num_coefs)
#     Q(Z) = log(Q0(Z,1,1,0,0))/(Z - S)
#     xp = as[1,1]*k*(-1.0+1.0im)
#     (S + k*as[1,1]) * exp(
#         (T(1.0)/(T(2)*pi*im)) * (
#             sum(Fun(Q, Segment(-maxZ, xp), num_coefs)) +
#             sum(Fun(Q, Segment(xp,-xp), num_coefs)) +
#             sum(Fun(Q, Segment(-xp,maxZ), num_coefs))
#         )
#     )
# end
