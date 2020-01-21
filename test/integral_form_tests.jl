
function test_integral_form()

    using EffectiveWaves

    # physical parameters
    θin = 0.0
    k=1.;
    ho = 2;

    medium = Medium(1.0,1.0+0.0im)
    ω = real(k*medium.c)
    # specie = Specie(ρ=0.1,r=1.0, c=0.5, volfrac=0.2)
    specie = Specie(ρ=0.1,r=1.0, c=0.5, volfrac=0.1)
    # specie = Specie(ρ=0.1,r=0.1, c=0.5, volfrac=0.1)
    specie2 = Specie(ρ=2.0, r=1.0, c=0.1, volfrac=0.15)

    # From effective wave theory
    k_eff0 = wavenumber_low_volumefraction(ω, medium, [specie]; basis_order = ho)
    k_effs = wavenumbers(ω, medium, [specie]; mesh_points = 10, tol = 1e-8, basis_order = ho)
    k_effs = sort(k_effs, by=imag)

    eff_medium = effective_medium(medium, [specie])
    ω/eff_medium.c
    k_effs2 = wavenumbers(ω, medium, [specie2]; tol = 1e-8, basis_order = ho)
    k_effs2 = sort(k_effs2, by=imag)

    using OffsetArrays, ApproxFun, IterTools
    include("src/integral_form/integral_form.jl")

    # (x, (MM_quad,b_mat)) = discrete_wave_system(ω, medium, specie; θin = θin, mesh_points = 501, basis_order=ho);
    X = 0.0:0.005:30.0
    (MM_quad,b_mat) = discrete_wave_system(ω, X, medium, specie;θin = θin, basis_order=ho);

    # discretization parameters
    J = length(collect(x)) - 1

    len = (J + 1) * (2ho + 1)
    MM_mat = reshape(MM_quad, (len, len));
    b = reshape(b_mat, (len));

    As = MM_mat\b;
    As_mat = reshape(As, (J+1, 2ho+1));

    amps_eff1 = scattering_amplitudes_average(ω, x, medium, [specie];
            k_eff = k_effs[1], basis_order = ho, θin=θin, tol=1e-8)
    amps_eff2 = scattering_amplitudes_average(ω, x, medium, [specie];
            k_eff = k_effs[2], basis_order = ho, θin=θin, tol=1e-8)
    amps_eff3 = scattering_amplitudes_average(ω, x, medium, [specie];
            k_eff = k_effs[end], basis_order = ho, θin=θin, tol=1e-8)

    amps0_eff = scattering_amplitudes_average(ω, x, medium, [specie];
            k_eff = k_eff0, basis_order = ho, θin=θin)
    amps2_eff = scattering_amplitudes_average(ω, x, medium, [specie2];
            k_eff = k_effs2[2], basis_order = ho, θin=θin, tol=1e-8)
    # sanity check: the abs of reflection coefficients should always be smaller than one.
    reflection_coefficient(ω, medium, specie2; amps = amps2_eff, θin = θin)
    reflection_coefficient(ω, medium, specie; amps = amps_eff2, θin = θin)
    reflection_coefficient(ω, medium, specie; amps = amps_eff3, θin = θin)
    reflection_coefficient(ω, medium, specie; amps = amps_eff1, θin = θin)

    error0_eff = reshape( abs.(MM_mat*amps0_eff.amplitudes[:] .- b), (J+1, 2ho+1))
    error_eff1 = reshape( abs.(MM_mat*amps_eff1.amplitudes[:] .- b), (J+1, 2ho+1))
    error_eff2 = reshape( abs.(MM_mat*amps_eff2.amplitudes[:] .- b), (J+1, 2ho+1))
    error_eff3 = reshape( abs.(MM_mat*amps_eff3.amplitudes[:] .- b), (J+1, 2ho+1))

    error2_eff = reshape( abs.(MM_mat*amps2_eff.amplitudes[:] .- b), (J+1, 2ho+1))
    # With a completely wrong wavenumber I can still get a small error by just scalling the amplitudes
    amps2_eff = scattering_amplitudes_average(ω, x, medium, [specie];
            k_eff = k_effs2[2], basis_order = ho, θin=θin)
    error2_eff = reshape( abs.(MM_mat*amps2_eff.amplitudes[:] .- b), (J+1, 2ho+1))
    # scale_amplitudes_effective(ω, k_effs[1], amps2_eff.amplitudes, medium, [specie]; θin = θin)

    error = reshape( abs.((MM_mat*As)./b .- 1.0+0.0im), (J+1, 2ho+1))

    using Plots; pyplot(linewidth=2) # unicodeplots()
    plot(xlabel = "depth (1 wavelength = 2π )", ylabel = "error %", ylims=(-0.1,0.5), title="Transmitted wave errors")
    plot!(x,error_eff1[:,ho+1], label = "Eff. error")
    plot!(x,error0_eff[:,ho+1], linestyle=:dash, label = "Eff. low φ error")
    plot!(x,error2_eff[:,ho+1], linestyle=:dot, label = "Eff. wrong k_eff error")
    plot!(x,error[:,ho+1], linestyle=:dashdot, label = "Integral method error")

    plot(x, [real.(As_mat[:,ho+1]),imag.(As_mat[:,ho+1])], labels = ["real sol." "imag sol."])
    plot!(x, [real.(amps_eff1.amplitudes[:,ho+1]),imag.(amps_eff1.amplitudes[:,ho+1])],
        labels = ["real eff. 1" "imag eff. 1"], linestyle=:dash)
    plot!(x, [real.(amps_eff2.amplitudes[:,ho+1]),imag.(amps_eff2.amplitudes[:,ho+1])],
        labels = ["real eff. 2" "imag eff. 2"], linestyle=:dash)

    plot(x, [abs.(As_mat[:,ho+1] .- amps_eff1.amplitudes[:,ho+1])],
            labels = ["diff abs eff. " "diff imag eff."])

    plot(x, [real.(As_mat[:,ho+2]),imag.(As_mat[:,ho+2])], labels = ["real sol." "imag sol."])
    plot!(x, [real.(amps_eff1.amplitudes[:,ho+2]),imag.(amps_eff1.amplitudes[:,ho+2])],
        labels = ["real eff. 1" "imag eff. 1"], linestyle=:dash)

    # plot!(x, [real.(amps_eff3.amplitudes[:,ho+1]),imag.(amps_eff3.amplitudes[:,ho+1])],
    #     labels = ["real eff. 3" "imag eff. 3"], linestyle=:dash)
    # plot!(x, [real.(amps0_eff.amplitudes[:,ho+1]),imag.(amps0_eff.amplitudes[:,ho+1])],
    #     labels = ["real φ eff." "imag φ eff."], linestyle=:dot)

    is = 250:(length(collect(x))-1)
    plot(x[is], log.(abs.(As_mat[is,ho+1])), labels = "abs sol.")
    plot!(x[is], log.(abs.(As_eff_mat[is,ho+1])), labels = "abs eff.")
    # b_mat ≈ [ sum(MM_quad[l,m,j,n]*A_mat[j,n] for j=1:(J+1), n=1:(2M+1)) for l=1:(J+1), m=1:(2M+1)]
end

# ints = [ check_integration(1.0+1.0im; h = 1./n) for n=10:30:310]

function test_check_integration()
#from Mathematics, had several errors when running
    math = [
            [65.0544 - 146.891im, 24.1498 - 8.37225im, -5.10012 + 3.74221im
                , 0.940428 - 3.4201im, 7.25312 - 0.0534466im, 22.9058 + 38.434im, -49.5574 + 92.4377im]
            ,[-32.7224 + 1.52201im, -12.033 - 16.2499im, 5.93884 - 10.8663im
                , 8.51838 + 3.68509im, -4.94285 + 9.32462im, -12.0994 - 5.48456im, 3.6755 - 13.9906im]
            ,[-9.82076 - 18.6418im, 10.2188 - 12.2782im, 11.8699 + 5.34705im
                , -2.698 + 12.0431im, -13.1158 - 0.218682im, -2.92199 - 13.8296im, 12.9993 - 6.4169im]
         ];
    # result of ints = check_integration(h = 1/400., max_x = 8.0)
    julia = [
            [64.9882-146.903im, 24.1561-8.37069im, -5.10105+3.74215im
                , 0.940588-3.4198im, 7.25286-0.0526824im, 22.9026+38.4402im, -49.5882+92.5011im]
            ,[-32.8359+1.66127im, -12.0424-16.2309im, 5.93782-10.8642im
                , 8.51816+3.68541im, -4.94155+9.32627im, -12.1151-5.49728im, 3.77515-13.8555im]
            ,[-9.88661-18.6492im, 10.2111-12.2772im, 11.8693+5.34716im
                , -2.69806+12.0429im, -13.1162-0.21808im, -2.92135-13.8369im, 12.9773-6.359im]
        ];

    using Plots; pyplot()

    scatter(-3:3, abs.(julia[1]./math[1] .- 1), xlab = "basis_order", ylab = "rel. error", lab ="", title="Mathematica vs Julia for x = 0.0")
    scatter(-3:3, abs.(julia[2]./math[2] .- 1), xlab = "basis_order", ylab = "rel. error", lab ="", title="Mathematica vs Julia for x = 1.0")
    scatter(-3:3, abs.(julia[3]./math[3] .- 1), xlab = "basis_order", ylab = "rel. error", lab ="", title="Mathematica vs Julia for x = 2.0")

#from julia
# using JLD
# ints =  first(values(load("integrated_As.jld")))
# data = [ [ints[i][j][k] for i in eachindex(ints)] for j=1:3, k=1:7];
# # i=2;j=3;k=5; ints[i][j][k] == data[j,k][i]
# # x = 0. => j = 1
# M = 3;
#
# hs =  [ 1./n for n=10:30:310]
# j = 3
#  plot()
#  for k=1:7
#      plot!(hs[4:end], abs.(data[j,k][4:end]), label = "m = $(k - M - 1)", xlab = "mesh element size")
#      scatter!([hs[end]], [abs(math[j][k])])
#  end
#   gui()
#
#  plot()
#  for k=1:7
#      plot!(hs[4:end], 1 - abs.(data[j,k][4:end])./abs(data[j,k][end]), label = "m = $(k - M - 1)")
#  end
#   gui()
end

# tests the integration scheme
function check_integration(k_eff::Complex{Float64} = 1.0+1.0im; k=1.,a=1., h = a*k/55., max_x = 6.0)

    A(n,x) = im^Float64(n)*exp(-im*n*θin)*exp(im*x*k_eff)
    k=1.; a=1.;
    # physical parameters
    θin = 0.3

    # discretization parameters
    M = 3;
    J = Int(round(max_x/h)) # choose an even number for Integration schemes
    x = OffsetArray{Float64}(undef, 0:J)
    x[0:J] = (0:J)*h

    # trapezoidal
    σ =  OffsetArray(trapezoidal_scheme(collect(x); xn=max_x), 0:J)

    PQ_quad = intergrand_kernel(x, a*k; θin = θin, M = M)

    # Apply integration scheme
    for j=0:J
        PQ_quad[:,:,j+1,:] *= σ[j]
    end
    # PQ_quad = [σ[j]*PQ_quad[l+1,m+M+1,j+1,n+M+1] for l=0:J, m=-M:M, j=0:J, n=-M:M]

    len = (J + 1) * (2M + 1)
    PQ_mat = reshape(PQ_quad, (len, len))

    A_mat = [ A(m,x[l]) for l = 0:J, m = -M:M]
    As = reshape(A_mat, (len))
    As_integrated = reshape(PQ_mat*As, (J+1,2M+1))

    i1 = Int(round(1/h)) # for x1 = 1
    i2 = Int(round(2/h)) # for x1 = 2
    [As_integrated[1,:], As_integrated[1+i1,:], As_integrated[1+i2,:]]
end
