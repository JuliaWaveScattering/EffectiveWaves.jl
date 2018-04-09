const opt_methods = (:adaptive_de_rand_1_bin, :adaptive_de_rand_1_bin_radiuslimited, :separable_nes, :xnes, :de_rand_1_bin, :de_rand_2_bin, :de_rand_1_bin_radiuslimited, :de_rand_2_bin_radiuslimited, :random_search, :generating_set_search, :probabilistic_descent)

# opt_fields = [(:r,(0.01,2.0)),(:volfrac,(0.,0.14)),(:ρ,(0.0,100.))]
# fix_fields = [(:c,[1.0+0.0im,1.0+0.0im])]
# optf = f_slow1
# medium = Medium{Float64}(1.0,1.)
# ωs = 0.01:0.01:1.0
# num_species=2
# MaxTime=5
# method = :adaptive_de_rand_1_bin

function optimal_species(optf, medium, ωs;
                opt_fields = [(:r,(0.1,2.0)), (:volfrac,(0.01,0.1))],
                fix_fields = [(:ρ,[0.0,Inf]), (:c,[1.0+0.0im,1.0+0.0im])],
                num_species=2, MaxTime=100., method=:choose_method,
                opt_kws...
)

  function F(; ρ=[0.,Inf], r=[0.5,0.5], c=[1.0+0.0im,1.0+0.0im], volfrac=[0.1,0.1])
    N = num_species
    sps = [Specie(ρ[i], r[i], c[i]; volfrac=volfrac[i]) for i=1:N]
    kTs = wavenumber_low_volfrac(ωs, medium, sps)
    # kTs = [wavenumber_low_volfrac(ω, medium, sps) for ω in ωs]
    optf(sps, medium, ωs, kTs)
  end
  function optF(x)
    var_fields = map(1:length(opt_fields)) do n
      (opt_fields[n][1], x[(1+(n-1)*num_species):(n*num_species)])
    end
    F(; var_fields..., fix_fields...)
  end
  if method == :choose_method
    fits =[
      best_fitness(bboptimize(optF; MaxTime = min(100.0,MaxTime), Method = m,
        SearchRange =  repeat([o[2] for o in opt_fields], inner=[num_species])))
    for m in opt_methods]
    method = opt_methods[findmin(fits)[2]]
    println("For MaxTime=$MaxTime, the optimal method was = $method.")
  end

  res = bboptimize(optF; MaxTime = MaxTime, Method=method,
    SearchRange =  repeat([o[2] for o in opt_fields], inner=[num_species]))
  parameters = reshape(best_candidate(res), num_species, length(opt_fields))
  # include the fixed paramters
  parameters = hcat(parameters, hcat([f[2] for f in fix_fields]...))
  fields = [f[1] for f in (opt_fields..., fix_fields...)]
  [ Specie(Float64; zip(fields, parameters[n,:])...) for n =1:num_species]
end
