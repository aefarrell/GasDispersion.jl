@testset "Simple turbulent jet tests" begin
    # example scenario
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

    # test type inheritance
    @test isa(plume(ex, SimpleJet()), Plume)

    # missing model params
    @test_throws MissingException plume(ambient, SimpleJet())

    # horizontal jet
    j = plume(ex, SimpleJet(release_angle=0.0, k2=a, k3=b))
    @test j(x,y,z) ≈ c
    @test j(-x,y,z) ≈ 0.0

end
