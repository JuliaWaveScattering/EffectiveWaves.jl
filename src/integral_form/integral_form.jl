# using ApproxFun
# using OffsetArrays
# using EffectiveWaves

function integrate_B_full(n::Int,X, Y0; Y1 =1000000, θin = 0.0, num_coefs = 10000)
    K(Y) = cos(Y*sin(θin) + n*atan2(Y,X))*hankelh1(n,sqrt(X^2+Y^2))
    # approximate function with Chebyshev polynomial (to high precision) then integrate from Y0 to Y1
    return 2.0*(-1.0)^n*sum(Fun(K,Y0..Y1, num_coefs))
end

# Y0 = sqrt(k^a12^2 - X^2)
function integrate_B(n::Int,X, Y0; θin = 0.0, num_coefs = 10000)
    Y1 = max(2000.0*X, 4000.0) # note Y1 is non-dimensional!
    # assymptotically approximate the integral from Y1 to Inf (tested in integrate_hankels.nb)
    Binf = (1.0+1.0im)*exp(im*Y1*(1.0 - sin(θin)))*
        (1.0 + (-1.0)^n*exp(2.0im*Y1*sin(θin))*(1.0 - sin(θin)) + sin(θin))/(sqrt(pi*Y1)*cos(θin)^2)

    return Binf + integrate_B_full(n, X, Y0; Y1=Y1, θin=θin, num_coefs = num_coefs)
end

# for only whole-correction, this doesn't involve an integral
function integrate_S(n::Int,X; θin = 0.0)
    S = 2.0*(im^Float64(n))*exp(-im*n*θin)*exp(im*X*cos(θin))/cos(θin)
    if X<0 S = conj(S) end
    S
end

function test_integral_form()

    # physical parameters
    θin = 0.0
    k=1.;
    medium = Medium(1.0,1.0+0.0im)
    ω = real(k*medium.c)
    specie = Specie(ρ=0.1,r=0.5, c=0.1, volfrac=0.15)
    specie2 = Specie(ρ=0.0,r=1.0, c=0.0, volfrac=0.15)

    (x, (MM_quad,b_mat)) = integral_form(k, medium, specie; θin = θin, mesh_points = 501);

    # discretization parameters
    M = Int( (size(b_mat,2) - 1)/2 )
    J = length(collect(x)) - 1

    # From effective wave theory
    k_eff0 = wavenumber_low_volfrac(ω, medium, [specie])
    k_eff = wavenumber(ω, medium, [specie])
    k_eff2 = wavenumber(ω, medium, [specie2])

    As0_fun = transmission_scattering_coefficients_field(ω, k_eff0, medium, [specie];
            max_hankel_order=M, θin=θin)
    As_fun = transmission_scattering_coefficients_field(ω, k_eff, medium, [specie];
            max_hankel_order=M, θin=θin)
    As2_fun = transmission_scattering_coefficients_field(ω, k_eff2, medium, [specie2];
            max_hankel_order=M, θin=θin)

    len = (J + 1) * (2M + 1)
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    As = MM_mat\b;

    As_eff_mat = transpose(hcat(As_fun.(x)...))
    As_eff = reshape(As_eff_mat, (len))
    As0_eff_mat = transpose(hcat(As0_fun.(x)...))
    As0_eff = reshape(As_eff_mat, (len))

    As2_eff_mat = transpose(hcat(As2_fun.(x)...))
    As2_eff = reshape(As2_eff_mat, (len))

    error0_eff = reshape( abs.((MM_mat*As0_eff)./b .- 1.0+0.0im), (J+1, 2M+1))
    error_eff = reshape( abs.((MM_mat*As_eff)./b .- 1.0+0.0im), (J+1, 2M+1))
    error2_eff = reshape( abs.((MM_mat*As2_eff)./b .- 1.0+0.0im), (J+1, 2M+1))

    As = MM_mat\b
    error = reshape( abs.((MM_mat*As)./b .- 1.0+0.0im), (J+1, 2M+1))

    using Plots; pyplot(linewidth=2)
    plot(xlabel = "depth (1 wavelength = 2π )", ylabel = "error %", ylims=(-0.1,1.5), title="Transmitted wave errors")
    plot!(collect(x),error_eff[:,M+1], label = "Eff. error")
    plot!(collect(x),error0_eff[:,M+1], linestyle=:dash, label = "Eff. low φ error")
    plot!(collect(x),error2_eff[:,M+1], linestyle=:dot, label = "Eff. wrong k_eff error")
    plot!(collect(x),error[:,M+1], linestyle=:dashdot, label = "Integral method error")

    As_mat = reshape(As, (J+1, 2M+1))

    plot(collect(x), [real.(As_mat[:,M+1]),imag.(As_mat[:,M+1])], labels = ["real sol." "imag sol."])
    plot!(collect(x), [real.(As_eff_mat[:,M+1]),imag.(As_eff_mat[:,M+1])], labels = ["real eff." "imag eff."])
    is = 250:(length(collect(x))-1)
    plot(collect(x[is]), log.(abs.(As_mat[is,M+1])), labels = "abs sol.")
    plot!(collect(x[is]), log.(abs.(As_eff_mat[is,M+1])), labels = "abs eff.")
    # b_mat ≈ [ sum(MM_quad[l,m,j,n]*A_mat[j,n] for j=1:(J+1), n=1:(2M+1)) for l=1:(J+1), m=1:(2M+1)]

    return x, As_mat
end

function integral_form_scattering_coefficients(ω::T,medium::Medium{T},specie::Specie{T};
        θin = 0.0)

    # physical parameters

    k=1.;
    medium = Medium(1.0,1.0+0.0im)
    k = ω/medium.c
    (x, (MM_quad,b_mat)) = integral_form(k, medium, specie; θin = θin, mesh_points = 501);

    # discretization parameters
    M = Int( (size(b_mat,2) - 1)/2 )
    J = length(collect(x)) - 1

    len = (J + 1) * (2M + 1)
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    As = MM_mat\b
    As_mat = reshape(As, (J+1, 2M+1))

    return x, As_mat
end

function reflection_coefficient_integrated(ω::T, medium::Medium, specie::Specie; A_mat::Array{Complex{T},2} = ) where T <: Number

As_mat = reshape(As, (J+1, 2M+1))

end

function integral_form(k::Float64, medium::Medium, specie::Specie;
        θin::Float64 = 0.0,
        x::AbstractVector = [0.], mesh_points::Int = 501,
        hankel_order = maximum_hankel_order(k*medium.c, medium, [specie]; tol=1e-3))

    ak = k*specie.r;
    ω = real(k*medium.c)
    M = hankel_order;

    # estimate a large enough mesh
    if x == [0.]
        k_eff = wavenumber_low_volfrac(ω, medium, [specie])
        max_x = 10.0*k/imag(k_eff) # at this A ≈ exp(-10) ≈ 4.5e-5
        J = mesh_points - 1
        h = ak/Int(round(J*ak/max_x));
        x = OffsetArray((0:J)*h, 0:J)
    else
        J = length(collect(x)) - 1
        h = x[2] - x[1]
    end

    Z = OffsetArray{Complex{Float64}}(-M:M);
    for m = 0:M
        Z[m] = Zn(ω,specie,medium,m)
        Z[-m] = Z[m]
    end

    # integration scheme: trapezoidal
    σ =  OffsetArray(trap_scheme(collect(x)), 0:J)
    PQ_quad = intergrand_kernel(x; ak = ak, θin = θin, M = M);

    MM_quad = [
        specie.num_density*Z[n]*σ[j]*PQ_quad[l+1,m+M+1,j+1,n+M+1] + k^2*( (m==n && j==l) ? 1.0+0.0im : 0.0+0.0im)
    for  l=0:J, m=-M:M, j=0:J, n=-M:M];

    b_mat = [ -k^2*exp(im*x[l]*cos(θin))*exp(im*m*(pi/2.0 - θin)) for l = 0:J, m = -M:M]

    return (x, (MM_quad,b_mat))
end

function intergrand_kernel(x::AbstractVector; ak::Float64 = 1.0, θin::Float64 = 0.0,
        M::Int = 2, num_coefs::Int = 10000)

    dx = x[2] - x[1]
    J = length(collect(x)) -1

    if !(typeof(x) <: OffsetArray)
        if J*dx != x[end] warn("Unexpected x = $x.") end
        x = OffsetArray((0:J)*dx, 0:J)
    end
    if !(Int(floor(ak/dx)) ≈ ak/dx)
        warn("There are no mesh points exactly on-top of the intergrands kinks. This could lead to poor accuracy.")
    end
    p = min(Int(floor(ak/dx)),J)
    X = OffsetArray((-J:J)*dx, -J:J)

    B = OffsetArray{Complex{Float64}}(-p:p, -2M:2M);
    for j = -p:p, m = -2M:2M
        if ak^2 -X[j]^2 < -dx^2 error("evaluating B in the wrong domain") end
        B[j,m] = integrate_B(m, X[j], sqrt(abs(ak^2 -X[j]^2)); θin = θin, num_coefs=num_coefs)
    end
    S = OffsetArray{Complex{Float64}}(-J:J, -2M:2M);
    for j = -J:J, m = -2M:2M
        S[j,m] = integrate_S(m, X[j]; θin = θin)
    end
    function intergrand(l,j,m,n)
        P = S[j-l,n-m]
        Q = (abs(j-l)<= p) ? (B[j-l,n-m] - S[j-l,n-m]) : 0.0+0.0im
        P + Q
    end

    intergrand_quad = [intergrand(l,j,m,n) for  l=0:J, m=-M:M, j=0:J, n=-M:M]

    return intergrand_quad
end


# plot(collect(x),[real(A_mat[:,M+1]),imag(A_mat[:,M+1])])
# plot(collect(x),abs.(A_mat[:,M+1]))
