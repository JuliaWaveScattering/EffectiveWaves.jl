using EffectiveWaves, Test, LinearAlgebra

# # Had to remove from testset block do to strange scoping error
# function randlm(order)
#     l = rand(0:order)
#     m = rand(-l:l)
#     return l, m
# end


@testset "fourth order matrix test" begin

    L = 2
    L1 = 3

    B = [
        ((l,m),(l2,m2),(dl,dm),(l1,m1))
                for dl = 0:L for dm = -dl:dl
            for l1 = 0:L1 for m1 = -l1:l1
        for l = 0:L for m = -l:l
    for l2 = 0:L1 for m2 = -l2:l2]

    len = (L1+1)^2 * (L+1)^2
    B = reshape(B, (:,len))

    eigF = [
        ((l,m),(l2,m2))
        for l = 0:L for m = -l:l
    for l2 = 0:L1 for m2 = -l2:l2]

    function randlm(order)
        local l = rand(0:order) # need unique names due to scoping
        local m = rand(-l:l)
        return l, m
    end

    l,m = randlm(L)
    l2,m2 = randlm(L1)

    dl,dm = randlm(L)
    l1,m1 = randlm(L1)

    n = lm_to_spherical_harmonic_index(l,m)
    n2 = lm_to_spherical_harmonic_index(l2,m2)
    dn = lm_to_spherical_harmonic_index(dl,dm)
    n1 = lm_to_spherical_harmonic_index(l1,m1)

    @test eigF == [b[3:4] for b in B[(n - 1) * (L1+1)^2 + n2,:]]
    @test eigF == [b[1:2] for b in B[:,(dn - 1) * (L1+1)^2 + n1]]

    @test ((l,m),(l2,m2),(dl,dm),(l1,m1)) == B[(n - 1) * (L1+1)^2 + n2, (dn - 1) * (L1+1)^2 + n1]

end
# @testset "MyTest" begin
#
#     function myfun(m)
#         n = m
#         return n
#     end
#
#     n = myfun(1)
#
#     println("n = ", n)
#
#     myfun(3) # <-- changes n to 3
#
#     println("n = ", n)
#
# end
