using ApproxFun

function integrate_B_full(n::Int,X, Y0; Y1 =1000000, θin = 0.0)
    K(Y) = cos(Y*sin(θin) + n*atan2(Y,X))*hankelh1(n,sqrt(X^2+Y^2))
    # approximate function with Chebyshev polynomial (to high precision) then integrate from Y0 to Y1
    return 2.0*(-1.0)^n*sum(Fun(K,Y0..Y1))
end

# assuming for now whole-correction
function S(n::Int,X; θin = 0.0)
    2.0*(im)^n*exp(-im*n*θin)*exp(im*X*cos(θin))/cos(θin)
end

# Y0 = sqrt(k^a12^2 - X^2)
function integrate_B(n::Int,X, Y0; θin = 0.0)
    Y1 = max(2000.0*X, 4000.0) # note Y1 is non-dimensional!
    # assymptotically approximate the integral from Y1 to Inf (tested in integrate_hankels.nb)
    Binf = (1.0+1.0im)*exp(im*Y1*(1.0 - sin(θin)))*
        (1.0 + (-1.0)^n*exp(2.0im*Y1*sin(θin))*(1.0 - sin(θin)) + sin(θin))/(sqrt(pi*Y1)*cos(θin)^2)

    return Binf + integrate_B_full(n, X, Y0; Y1=Y1, θin=θin)
end

function integral_form(ks)

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
    sum(S_fun*test_fun) ≈ 2*(exp(10.*(im*cos(θin)-1))-1.0)/(cos(θin)*(im*cos(θin)-1.0))

    Sop = DefiniteIntegral(Chebyshev(d))[S_fun]
    Sop*test_fun ≈ 2*(exp(10.*(im*cos(θin)-1))-1.0)/(cos(θin)*(im*cos(θin)-1.0))

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
    BS_fun = Fun((x,X) -> (abs(X-x) < 1.0) ? B_fun(X-x) - S_fun(X-x) : 0.0+0.0im, d^2,100000)

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
    BS_op = DefiniteIntegral(Legendre(d))[
        LowRankFun((x,X) -> (abs(X-x) < 1.0) ? cos(pi*(X-x)/2.0): 0.0+0.0im, d^2)
    ]
    BS_op = DefiniteIntegral(Legendre(d))[
        Fun((x,X) -> cos(pi*(X-x)/2.0), d^2,100)
    ]


    one = Fun(x->1.0, d)
    u = (I - V)\one

    # Test
    x = Fun(d)
    u(0.1)-cumsum(exp(-(0.1-x))*u)(2.0) # ≈ 1

end
