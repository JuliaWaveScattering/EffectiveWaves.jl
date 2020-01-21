using ApproxFun
using Plots; pyplot()

include("numerical_integration.jl")


x1=0.5; x0=1.5; x2=2.;α=4.;

f(x) = (x<x1) ? x - x1 - exp(-α*(x1-x0)^8) : ( (x<x2) ? -exp(-α*(x-x0)^2.0) : -exp(-α*(x2-x0)^2.0) -(x-x2)  )
A(x) = cos(1.5*x)
g(x) = f(x)*A(x)

y1 = x1-1.0
y2 = x2+1.0
plot(g,y1:0.01:y2)

g_fun = Fun(g,y1..y2)
sumg = sum(g_fun)

plot!(g_fun,y1:0.01:y2)

# ns = 4:4:20
ns = 100:10:190
ns = [ns; 200:100:1000]
ns = [ns; 2000:1000:10000]

# puts points on the discontinuities
xs = map(ns) do n
    nxs = Int(round(n*(x2 - x1)/(y2-y1)))
    if isodd(nxs) nxs = nxs + 1 end
    h = (x2 - x1)/nxs
    x = reverse(collect(x1:-h:y1))
    x = [x; collect(x1+h:h:y2)]
end

hs = [x[2] - x[1] for x in xs]
plot(hs, [sumg])

ints_trap = map(hs) do h
    x = y1:h:y2
    σ = trapezoidal_scheme(x; x0 = y1, xn = y2)
    sum(g(x[i])*σ[i] for i in eachindex(x))
end

ints_simp = map(xs) do x
    σ = simpson_scheme(x; x0 = y1, xn = y2)
    sum(g(x[i])*σ[i] for i in eachindex(x))
end


scatter!(hs, ints_trap)
plot!(hs, ints_trap)

scatter!(hs, ints_simp)
plot!(hs, ints_simp)
