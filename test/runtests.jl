using SNMRForward, Test, LinearAlgebra

const zero_tol = 1e-15

@testset "Magnetic field decompositions" begin
    # perpendicular magnetic fields
    @test isapprox(SNMRForward.perpendicular_field(1., 1., π/2, 0.), [1., 0.])
    @test isapprox(SNMRForward.perpendicular_field(1., 1., 0., 0.), [-1., 0.])
    @test norm(SNMRForward.perpendicular_field(1., 0., π/2, 0.)) < zero_tol

end