using ApproxFun
using HCubature


a12 =1.0

n=1
k=1.
X = 0.2

function Bfull(n::Int,X; Y0= sqrt(k^2*a12^2-X^2), Y1 =1000000 ,θin = 0.0)
    K(Y) = cos(Y*sin(θin) + n*atan2(Y,X))*hankelh1(n,sqrt(X^2+Y^2))
    return 2.0*(-1.0)^n*sum(Fun(K,Y0..Y1))
end

function B(n::Int,X; Y0= sqrt(k^2*a12^2-X^2), θin = 0.0)
    Y1 = max(2000.0*X, 4000.0)
    # assymptotically approximate the integral from Y1 to Inf
    Binf = (1.0+1.0im)*exp(im*Y1*(1.0 - sin(θin)))*
        (1.0 + (-1.0)^n*exp(2.0im*Y1*sin(θin))*(1.0 - sin(θin)) + sin(θin))/(sqrt(pi*Y1)*cos(θin)^2)

    return Binf + Bfull(n, X ;Y0=Y0, Y1=Y1, θin=θin)
end

K(Y) = cos(Y*sin(θin) + n*atan2(Y,X))*hankelh1(n,sqrt(X^2+Y^2))

interval = sqrt(k^2*a12^2-X^2)..1000000
interval1 = interval.left..1000.5
interval2 = interval1.right..interval.right
Ys = Fun(identity,interval)
Y1s = Fun(identity,interval1)
Y2s = Fun(identity,interval2)
Ks = K.(Ys)


sum(Ks) - sum(K.(Y1s)) - sum(K.(Y2s))
@time 2.0*(-1.0)^n*sum(Ks)
@time 2.0*(-1.0)^n.*hcubature(Y -> K(Y[1]), (interval.left,), (interval.right,))


Σ = DefiniteIntegral(Chebyshev(interval))
@time sum(Σ[i]*Ks.coefficients[i] for i in eachindex(Ks.coefficients))


interval = 0.0..10.0
α=1.0
x = Fun(interval)
Σ = DefiniteIntegral(Chebyshev(interval))

b = exp((im-α)*x)-1.0im/(im-α)
u = (I+1.0im*Σ)\b;

#Analytically the solution should be
norm(u.(0.0:0.1:5.0)-exp((im-α).*(0.0:0.1:5.0)))/norm(exp((im-α).*(0.0:0.1:5.0)))

@time sum(Σ[i]*Ks.coefficients[i] for i in eachindex(Ks.coefficients))


2*(-1)^n
