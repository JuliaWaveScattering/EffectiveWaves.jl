using ApproxFun, Plots
using LinearAlgebra, SpecialFunctions
# using HCubature

θin = 0.4; X = 1.0; n = 3
num_coefs = 1000;
Y0 = 0.0; Y1 = 2.0;

d = Domain(0..2)

K2D(X,Y) = cos((Y+X*im)*sin(θin) * n*dot(Y,X))*sin(sqrt(X^2+Y^2+ X*Y*im))
# K(Y) = cos((Y+X*im)*sin(θin) * n*dot(Y,X))*hankelh1(n,sqrt(X^2+Y^2+ X*Y*im))
K2Dreal(X,Y) = real(cos((Y+X*im)*sin(θin) * n*dot(Y,X))*sin(sqrt(X^2+Y^2+ X*Y*im)))
K2Dimag(X,Y) = imag(cos((Y+X*im)*sin(θin) * n*dot(Y,X))*sin(sqrt(X^2+Y^2+ X*Y*im)))

freal = Fun(K2Dreal, d^2, num_coefs,num_coefs)
fimag = ProductFun(K2Dimag, d^2, num_coefs,num_coefs)

f = Fun(K, Domain(0..2), num_coefs)
freal = Fun(Kreal,Domain(0..2), num_coefs)
fimag = Fun(Kreal,Domain(0..2), num_coefs)


function Bfull(n::Int,X; Y0= sqrt(k^2*a12^2-X^2), Y1 =1000000 , θin = 0.0)
    K(Y) = cos(Y*sin(θin) + n*atan(Y,X))*hankelh1(n,sqrt(X^2+Y^2))
    return 2.0*(-1.0)^n*sum(Fun(K,Y0..Y1))
end

function B(n::Int,X; Y0= sqrt(k^2*a12^2-X^2), θin = 0.0)
    Y1 = max(2000.0*X, 4000.0)
    # assymptotically approximate the integral from Y1 to Inf
    Binf = (1.0+1.0im)*exp(im*Y1*(1.0 - sin(θin)))*
        (1.0 + (-1.0)^n*exp(2.0im*Y1*sin(θin))*(1.0 - sin(θin)) + sin(θin))/(sqrt(pi*Y1)*cos(θin)^2)

    return Binf + Bfull(n, X ;Y0=Y0, Y1=Y1, θin=θin)
end


θin = 0.1
a12 =1.0
n=1
k=1.
X = 0.2
K(Y) = cos(Y*sin(θin) + n*atan(Y,X))*hankelh1(n,sqrt(X^2+Y^2))

interval = sqrt(k^2*a12^2-X^2)..1000000
interval1 = interval.left..1000.7
interval2 = interval1.right..interval.right
Ys = Fun(identity,interval)
Y1s = Fun(identity,interval1)
Y2s = Fun(identity,interval2)
Ks = K.(Ys)

interval1 = 0.0..1.0
B_app =  Fun(x -> B(0,x),interval1)

sum(Ks) - sum(K.(Y1s)) - sum(K.(Y2s))
@time 2.0*(-1.0)^n*sum(Ks)
@time 2.0*(-1.0)^n.*hcubature(Y -> K(Y[1]), (interval.left,), (interval.right,))

Σ = DefiniteIntegral(Chebyshev(interval))
@time sum(Σ[i]*Ks.coefficients[i] for i in eachindex(Ks.coefficients))


# Solving fredholm https://github.com/JuliaApproximation/ApproxFun.jl/issues/570
d = Interval(0.0, 2.0)
K2 = Fun((x,y) -> exp(-(x-y)), d^2)
V = DefiniteIntegral(Chebyshev(d))[LowRankFun(K2,d^2)]
one = Fun(x->1.0, d)
u = (I - V)\one

# Test
x = Fun(d)
u(0.1)-cumsum(exp(-(0.1-x))*u)(2.0) # ≈ 1

# Solving simple fredholm http://juliaapproximation.github.io/ApproxFun.jl/latest/usage/equations.html#Linear-equations-1
interval = 0.0..10.0
α=1.0
x = Fun(interval)
Σ = DefiniteIntegral(Chebyshev(interval))

b = exp((im-α)*x)-1.0im/(im-α)
u = (I+1.0im*Σ)\b;

#Analytically the solution should be
norm(u.(0.0:0.1:5.0)-exp((im-α).*(0.0:0.1:5.0)))/norm(exp((im-α).*(0.0:0.1:5.0)))

@time sum(Σ[i]*Ks.coefficients[i] for i in eachindex(Ks.coefficients))

a =1.0
θin = 0.1
n=1
k=1.
X = 0.2

θin = 0.3

d = Interval(0.0, 10.0)
xs = Fun(d)
test_fun = exp(-xs)
S_fun = S(0,xs; θin = 0.3)
sum(S_fun*test_fun) ≈ 2*(exp(10 .* (im*cos(θin)-1))-1.0)/(cos(θin)*(im*cos(θin)-1.0))

Sop = DefiniteIntegral(Chebyshev(d))[S_fun]
Sop*test_fun ≈ 2*(exp(10 .* (im*cos(θin)-1))-1.0)/(cos(θin)*(im*cos(θin)-1.0))

x1 = Fun(Interval(0.0, 1.0))
@time f1 = Fun( x-> exp(-x), 0.0..1.0)
@time f3 = Fun( x-> (x>1.0) ? 0.0 : exp(-x), 0.0..10.0)
f2 = exp(-xs)

f4 = f1+f2

f1(0.3)+f2(0.3)
2*f4(0.3)

f1(1.3)+f2(1.3)
f4(1.3)

dd_fun = Fun((x,X) -> (abs(X-x) < 1.0) ? S_fun(0,X-x) :0.0+0.0im, d^2,10)

dd_fun = Fun((x,X) -> S_fun(X), d^2,10)

d1 = 0.0..1.0
B_fun = Fun(x -> integrate_B(0,x, sqrt(1.0 - x^2)), d1,40)#450)
S_fun = S(0,xs)

BS_lowfun = LowRankFun((x,X) -> (abs(X-x) < 1.0) ? B_fun(X-x) - S_fun(X-x) : 0.0+0.0im, d^2)
BS_lowfun = LowRankFun((x,X) -> (abs(X-x) < 1.0) ? B_fun(X-x) - S_fun(X-x) : 0.0+0.0im, d^2)
BS_fun = Fun((x,X) -> (abs(X-x) < 1.0) ? B_fun(X-x) - S_fun(X-x) : 0.0+0.0im, d^2,1000)

# r1 = 10.0*rand(10);
# r2 = 10.0*rand(10);
#
# m = mean( abs((abs(X-x) < 1.0) ? B_fun(X-x) - S_fun(X-x) : 0.0+0.0im)  for x in r1, X in r2)
# [ abs( BS_fun(x,X) - ((abs(X-x) < 1.0) ? B_fun(X-x) - S_fun(X-x) : 0.0+0.0im))  for x in r1, X in r2]
# [BS_fun(x,X) for x in r1, X in r2]

# Y0 = sqrt(k^a12^2 - (X-x)^2)
BS_op = DefiniteIntegral(Legendre(d))[
    BS_fun
]
# Y0 = sqrt(k^a12^2 - (X-x)^2)
BS_op = DefiniteIntegral(Chebyshev(d))[
    LowRankFun((x,X) -> (abs(X-x) < 1.0) ? B_fun(X-x) - S_fun(X-x) : 0.0+0.0im, d^2)
]

d = Interval(0.0, 2.0)
fun = Fun((x,X) -> (abs(X-x) < 1.0) ? exp(-abs(X-x)) : 0.0+0.0im, d^2,1000)
DefiniteIntegral(Chebyshev(d))[fun]
DefiniteIntegral(Legendre(d))[Fun((x,X) -> (abs(X-x) < 1.0) ? exp(-abs(X-x)) : 0.0+0.0im, d^2,1000)]
DefiniteIntegral(Chebyshev(d))[Fun((x,X) -> (abs(X-x) < 1.0) ? exp(-abs(X-x)) : 0.0+0.0im, d^2,1000)]
DefiniteIntegral(Chebyshev(d))[LowRankFun((x,X) -> (abs(X-x) < 1.0) ? exp(-abs(X-x)) : 0.0+0.0im, d^2)]
DefiniteIntegral(Chebyshev(d))[LowRankFun((x,X) -> exp(-abs(X-x)), d^2)]


f = LowRankFun((x,X) -> exp(-abs(X-x)))
f1 = Fun(x -> exp(-abs(x)),d,100)
f = Fun((x,X) -> exp(-abs(X-x)),d^2,100)

int_op = DefiniteIntegral(Legendre(d))[
LowRankFun((x,X) -> (abs(X-x) < 1.0) ? cos(pi*(X-x)/2.0): 0.0+0.0im, d^2)
]
BS_op = DefiniteIntegral(Legendre(d))[
Fun((x,X) -> cos(pi*(X-x)/2.0), d^2,100)
]
