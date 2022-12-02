

"""
    asymptotic_monopole_wavenumbers(ω, medium::Acoustic, species; num_wavenumbers = 2)

Calculates the asymptotic effective wavenumbers for monopole scatterers by assuming a large number of wavenumbers. The 2D results are deduced in [Section 5](https://arxiv.org/pdf/1905.06996.pdf).
"""
asymptotic_monopole_wavenumbers(ω::T, medium::Acoustic{T}, sps::Species; kws...) where T = asymptotic_monopole_wavenumbers(ω, Microstructure(medium,sps); kws...)

function asymptotic_monopole_wavenumbers(ω::T, micro::Microstructure{2};
        num_wavenumbers = 2
    ) where T

    medium = micro.medium
    species = micro.species

    P = Int(round(num_wavenumbers / 2) + 1)
    k = ω / medium.c

    b = mean(s1.separation_ratio * outer_radius(s1) + s2.separation_ratio * outer_radius(s2) for s1 in species, s2 in species)
    num_density = mean(number_density(s) for s in species)

    T0 = mean(t_matrix(s.particle, medium, ω, 0) for s in species)[1]

    c = sqrt(2pi) * num_density * b^2 * T0 * hankelh1(0,k * b) * exp(-im * pi / 4.0)

    rc = abs(c)
    θc = angle(c)

    kps = map(-ceil(θc / (2pi)):P) do p
        σ = θc + 2pi * p
        kp = (σ + im * log(abs(σ)^(3/2) / rc)) / b
    end

    kns = map(ceil(θc / (2pi) - 3 / 4):P) do p
        σ = θc - 3pi / 2 - 2pi * p
        kn = (σ + im * log(abs(σ)^(3/2) / rc)) / b
    end

   return [kns; kps]
end

function asymptotic_monopole_wavenumbers(ω::T, micro::Microstructure{3};
        num_wavenumbers = 2
    ) where T

    medium = micro.medium
    species = micro.species

    P = Int(round(num_wavenumbers / 2) + 1)
    k = ω / medium.c

    b = mean(s1.separation_ratio * outer_radius(s1) + s2.separation_ratio * outer_radius(s2) for s1 in species, s2 in species)
    num_density = mean(number_density(s) for s in species)

    T0 = mean(t_matrix(s.particle, medium, ω, 0) for s in species)[1]

    c = 2pi * num_density * b^2 * T0 * shankelh1(0,k * b)

    rc = abs(c)
    θc = angle(c)

    kps = map(-P:P) do p
        σ = θc + 2pi * p
        kp = (σ + im * log(σ ^ 2 / rc)) / b
    end

   return kps
end
