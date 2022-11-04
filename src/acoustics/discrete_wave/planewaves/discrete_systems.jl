"note that this uses the non-dimensional X = k*depth"
function discrete_wave_system(ω::T, X::AbstractVector{T}, source::PlaneSource{T,2,1,Acoustic{T,2}}, material::Material{2,Halfspace{T,2}};
        tol::T = 1e-6,
        scheme::Symbol = :trapezoidal,
        basis_order::Int = 2,
        kws...
    ) where T<:AbstractFloat

    specie = material.microstructure.species[1]
    t_vec = t_matrix(specie, source.medium, ω, basis_order)

    k = real(ω / source.medium.c)
    a12k = specie.separation_ratio * T(2)*real(k * outer_radius(specie));
    M = basis_order;

    θin = transmission_angle(source,material)

    J = length(X) - 1
    h = X[2] - X[1]

    PQ_quad = intergrand_kernel(X, a12k; M = M, θin = θin, scheme=scheme);

    MM_quad = [
        (number_density(specie) / (k^2))*t_vec[m+M+1,m+M+1]*PQ_quad[l,m+M+1,j,n+M+1] - ( (m==n && j==l) ? 1.0+0.0im : 0.0+0.0im)
    for  l=1:(J+1), m=-M:M, j=1:(J+1), n=-M:M];

    b_mat = [
        -t_vec[m+M+1,m+M+1]*exp(im*X[l]*cos(θin))*exp(im*m*(pi/2.0 - θin))
    for l = 1:(J+1), m = -M:M]

    return (MM_quad,b_mat)
end

"Returns x the mesh used to discretise the integral equations."
function x_mesh(wave_eff_long::EffectivePlaneWaveMode{T}, wave_eff_short::EffectivePlaneWaveMode{T} = wave_eff_long;
        tol::T = T(1e-5),  a12::T = T(Inf),
        max_size::Int = 1000,
        min_size::Int = 5,
        max_x::T = -log(tol) / abs(imag(wave_eff_long.wavenumber))
        # max_x::T = -log(tol) / abs(imag(wave_eff_short.wavenumber))
        # max_x::T = (-log(tol))/abs(cos(wave_eff_long.θ_eff)*imag(wave_eff_long.wavenumber))
) where T<:AbstractFloat

    #= The default max_x results in:
        abs(exp(im*max_x*cos(θ_effs[end])*k_effs[end]/k)) < tol
    =#
    # estimate a reasonable derivative based on more rapidly varying wave_eff_short.
    df = abs(wave_eff_short.wavenumber)
    # NOTE was: df = abs(wave_eff_short.wavenumber * cos(wave_eff_short.θ_eff))

    # Based on Simpson's rule
        # dX  = (tol*90 / (df^4))^(1/5)
    # Based on trapezoidal integration
        dx  = (tol * T(30) / (df^2))^(1/3)

    # prioritise a mesh fine enough to give accurate discrete integrals
    if dx > a12 dx = a12 end

    if max_x/dx + 1 > max_size
        if dx == a12
            max_x = (max_size - 1)*dx
            @warn("The mesh max_size = $max_size which was too small for tol = $tol. Will shrink meshed region.")
        else # otherwise priortise shrink max_x and increase dx proportionally
            a = sqrt(max_x/((max_size - 1)*dx))
            dx = min(dx*a,a12)
            max_x = (max_size - 1)*dx
            @warn("The mesh max_size = $max_size which was too small for tol = $tol. Will make a smaller and coarser mesh.")
        end
    elseif max_x/dx + 1 < min_size
        max_x = (min_size - 1)*dx
    end

    # if whole correction length a12k was given, then make dX/a12k = integer
    if a12 != T(Inf) && dx < a12
        n = ceil(a12 / dx)
        dx = a12/n
    end

    return zero(T):dx:max_x
end
