using GSL
using Test, LinearAlgebra

@testset "Special functions" begin

    @testset "complex radial coordiantes" begin

    # Test 3-dimensional transforms
        xs = [rand(-1.01:0.1:1.0,3) + rand(-1.01:0.1:1.0,3)*im for i = 1:100]
        rθφs = cartesian_to_radial_coordiantes.(xs)
        @test maximum(norm.(xs - radial_to_cartesian_coordiantes.(rθφs))) < 2e-14

        xs = [rand(-1.01:0.1:1.0,3) for i = 1:100]
        rθφs = cartesian_to_radial_coordiantes.(xs)

        @test maximum(rθφ[2] for rθφ in rθφs) <= pi
        @test minimum(rθφ[2] for rθφ in rθφs) >= 0.0

        @test pi/2 < maximum(rθφ[3] for rθφ in rθφs) <= pi
        @test -pi <= minimum(rθφ[3] for rθφ in rθφs) < -pi/2

    # Test 2-dimensional transforms
        xs = [rand(-1.01:0.1:1.0,2) + rand(-1.01:0.1:1.0,2)*im for i = 1:100]
        rθs = cartesian_to_radial_coordiantes.(xs)

        @test maximum(norm.(xs - radial_to_cartesian_coordiantes.(rθs))) < 1e-14
    end
end
