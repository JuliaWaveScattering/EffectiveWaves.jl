# """
#     asymptotic_monopole_wavenumbers(ω, medium::PhysicalMedium, species; num_wavenumbers = 2)

# Estimates a box in complex k_eff space where all needed effective wavenumbers, for the given tolerance, are inside.
# """ # Should really use assymptotics for monopole scatterers to properly estimate box size. Would be more robust.
# function asymptotic_monopole_wavenumbers(ω::T, medium::PhysicalMedium, species::Species;
#         num_wavenumbers = 2
#     ) where T

#     # eff_medium = effective_medium(medium, species)
#     # ko = real(ω / eff_medium.c)
#     # if isnan(ko) ko = real(ω / medium.c) end

#     ko = real(ω / medium.c)

#     # After two wavelengths distance the solution would be smaller than tol
#     max_imag = - log(tol) * ko / (T(4) * pi)
#     max_real = T(4) * ko

#    return [[-max_real,max_real],[min_imag, max_imag]]
# end


"""
    box_keff(ω::T, medium::PhysicalMedium, species::Species;
            tol::T = 1e-5, min_imag::T = zero(T))

Estimates a box in complex k_eff space where all needed effective wavenumbers, for the given tolerance, are inside.
""" # Should really use assymptotics for monopole scatterers to properly estimate box size. Would be more robust.
function box_keff(ω::T, medium::PhysicalMedium, species::Species;
        tol::T = 1e-5, min_imag::T = zero(T)
    ) where T

    # eff_medium = effective_medium(medium, species)
    # ko = real(ω / eff_medium.c)
    # if isnan(ko) ko = real(ω / medium.c) end

    ko = real(ω / medium.c)

    # After two wavelengths distance the solution would be smaller than tol
    max_imag = - log(tol) * ko / (T(4) * pi)
    max_real = T(4) * ko

   return [[-max_real,max_real],[min_imag, max_imag]]
end

struct MySimplexer{T<:AbstractFloat} <: Optim.Simplexer
    dx::T
    dy::T
end

function simplexer(A::MySimplexer, initial_x::Array{T}) where T
    initial_simplex = [initial_x, initial_x + [zero(T), A.dy], initial_x + [A.dx, zero(T)]]
end

NelderMeadparameters(;α = 1.0, β = 1.0, γ = 1.0, δ = 0.5)::Optim.FixedParameters = Optim.FixedParameters(
    α = α, # reflection
    β = β, # expansion 1.0 = does not expand
    γ = γ, # contraction
    δ = δ  # shrink step
)

reduce_kvecs(vec::Vector,tol) = vec

function reduce_kvecs(vecs::Vector{Vector{T}},tol::T) where T<:AbstractFloat
    # All wavenumbers should attenuate
    inds = findall([w[2] < -tol for w in vecs])
    vecs[inds] = - vecs[inds]

    all_inds = collect(eachindex(vecs))
    vecs = map(vecs) do w
        ind_ins = findall([norm(v - w) < tol for v in vecs[all_inds]])
        inds = all_inds[ind_ins]
        deleteat!(all_inds,ind_ins)
        isempty(inds) ? [zero(T),-one(T)] :  mean(vecs[inds])
    end
    vecs = deleteat!(vecs, findall(v-> v == [zero(T),-one(T)], vecs) )

    return vecs
end

function reduce_kvecs(ks::Vector{Complex{T}},tol::T) where T<:AbstractFloat
    # All wavenumbers should attenuate
    inds = findall(imag.(ks) .< -tol)
    ks[inds] = - ks[inds]

    all_inds = collect(eachindex(ks))
    ks = map(ks) do w
        ind_ins = findall([norm(v - w) < tol for v in ks[all_inds]])
        inds = all_inds[ind_ins]
        deleteat!(all_inds,ind_ins)
        isempty(inds) ? zero(T)-one(T)*im :  mean(ks[inds])
    end

    ks = deleteat!(ks, findall(ks .== zero(T)-one(T)*im))

    return ks
end
