@testset "Simple turbulent jet tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 97
    ex = Scenario( Dict([
        :mass_emission_rate => 1.0,
        :jet_diameter => 1.0,
        :jet_velocity => 1.0,
        :jet_density => 1.2,
        :ambient_density => 1.2,
        :release_height => 1.0,
    ]))
    # known answers
    a, b = 6.0, 5.0
    ξ² = log(2)/b^2
    x, y, z =  a, √(a^2*ξ²), 1.0
    c = (2/π)*(1+exp(-25/9))

    # missing model params
    @test_throws MissingException plume(ambient, model=:simplejet)

    # horizontal jet
    j1 = plume(ex, model=:simplejet, release_angle=0.0, k₁=a, k₂=b)
    @test j1(x,y,z) ≈ c

end
