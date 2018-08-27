n = 8
function complex_least(n)
    w = rand(Complex{T},n)
    b = rand(Complex{T})
    U = rand(Complex{T},n,n)
    M = ctranspose(U)*U

    λ = b/(ctranspose(w)*inv(M)*w)
    α = λ*conj(inv(M))*conj(w)

    return (transpose(w)*α - b, transpose(α)*M*conj(α),α)
end

function complex_least_sys(i::Int)
    n=80; m=80;
    w = rand(Complex{T},m)
    c = rand(Complex{T},n)
    b = rand(Complex{T})
    A = rand(Complex{T},n,m)
    invM = inv(ctranspose(A)*A)

    λ = (b - transpose(w)*invM*ctranspose(A)*c)/(transpose(w)*invM*conj(w))

    α = λ*invM*conj(w) + invM*ctranspose(A)*c
    # α = invM*ctranspose(A)*c

    norm(transpose(w)*α - b)/norm(b)

    return norm(A*α - c)/norm(c)
end

# the wrong solution
function g(i::Int)
    n=80; m=80;
    w = rand(Complex{T},m)
    c = rand(Complex{T},n)
    b = rand(Complex{T})
    A = rand(Complex{T},n,m)
    invM = inv(ctranspose(A)*A)

    λ = (b - transpose(w)*invM*ctranspose(A)*c)/(transpose(w)*invM*(w))

    α = λ*invM*(w) + invM*ctranspose(A)*c
    # α = invM*ctranspose(A)*c

    norm(transpose(w)*α - b)/norm(b)
    norm(A*α - c)/norm(c)
end
mean(complex_least_sys(i) for i=1:500)
mean(g(i) for i=1:500)

plot([real.(A*α),real.(c)])

return norm(A*α - c)/norm(c)

ns = 2:2:1050
fs = zeros(Float64,length(ns))

for i in eachindex(ns)
    (constr, f, a) = complex_least(ns[i])
    fs[i] = abs(f)
end
# fs[2:end] = (fs[2:end]+fs[1:(length(fs)-1)])/2
v = exp.(-(-2:0.02:2).^2)
v = v/sum(v)
fs2 = conv(fs,v)
ns2 = conv(collect(ns),v)
fs2 = fs2[length(v):(length(fs2)-length(v))]
ns2 = ns2[length(v):(length(ns2)-length(v))]
# plot([fs2, fs[4:(length(fs)-3)]])
# plot([ns2, ns[4:(length(ns)-3)]])

plot(ns2,log.(fs2))
