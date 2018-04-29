using ApproxFun

x1=0.5; x0=1.5; x2=2.;α=4.;

f(x) = (x<x1) ? x - x1 - exp(-α*(x1-x0)^8) : ( (x<x2) ? -exp(-α*(x-x0)^2.0) : -exp(-α*(x2-x0)^2.0) -(x-x2)  )
A(x) = cos(1.5*x)
g(x) = f(x)*A(x)

y1 = x1-1.0
y2 = x2+1.0
plot(g,y1:0.01:y2)

g_fun = Fun(g,y1..y2)
sumf = sum(g_fun)

plot!(g_fun,y1:0.01:y2)

l = x1 - y1
hs = [l/i for i = 40:10:100]
hs = [hs; [l/i for i = 100:100:1000]]
hs = [hs; [l/i for i = 1000:1000:10000] ]

ints = map(hs) do h
    x = y1:h:y2
    sum(g(x[i])*h for i in eachindex(x))
end

plot(hs, [sumf])
scatter!(hs, ints)
plot!(hs, ints)

ints_trap = map(hs) do h
    x = y1:h:y2
    -(g(x[1]) + g(x[end]))*h/2.0 + sum(g(x[i])*h for i in eachindex(x))
end

# scatter!(hs, abs.(sumf .- ints_trap))
# plot!(hs, abs.(sumf .- ints_trap))
scatter!(hs, ints_trap)
plot!(hs, ints_trap)


function integrate_scheme(x::Vector{Float64})
    σs = [isodd(j) ? 4.0 : 2.0 for j in eachindex(x)]
    σs[1] = σs[end] = 1.0
    σs = σs*(x[2]-x[1])/3.0
    σs
end


ints_simp = map(hs) do h
    x = y1:h:y2
    σ = integrate_scheme(x)
    sum(g(x[i])*σ[i] for i in eachindex(x))
end

scatter!(hs, ints_simp)
plot!(hs, ints_simp)






ints_simp_3 = map(hs) do h
    x = y1:h:y2
    i1 = findmin(abs.(x-x1))[2]
    i2 = findmin(abs.(x-x2))[2]
    if( x[i1] != x1 || x[i2] != x2) error("did not find discontinuity points for h=$h and x=$x") end
    σ1 = integrate_scheme(x[1:i1])
    σ2 = integrate_scheme(x[i1:i2])
    σ1[end] += σ2[1]
    deleteat!(σ2,1)
    σ3 = integrate_scheme(x[i2:end])
    σ2[end] += σ3[1]
    deleteat!(σ3,1)
    σ = [σ1; σ2; σ3]
    sum(f(x[i])*σ[i] for i in eachindex(x))
end

scatter!(hs,ints_simp_3)
plot!(hs,ints_simp_3)
