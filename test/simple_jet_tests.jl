@testset "Simple turbulent jet tests" begin
    # example scenario
    ex = Scenario( Substance(name=:test,gas_density=1.2,liquid_density=1000.0,
                    boiling_temp=NaN,latent_heat=NaN,gas_heat_capacity=NaN,
                    liquid_heat_capacity=NaN),
                   Release(mass_rate = 1.0, duration = Inf, diameter = 1.0,
                            velocity = 1.0, height = 1.0, pressure = 0,
                            temperature = 0, fraction_liquid=0.0),
                   Ambient(density=1.2))
    # known answers
    a, b = 6.0, 5.0
    ξ² = log(2)/b^2
    x, y, z =  a, √(a^2*ξ²), 1.0
    c = (2/π)*(1+exp(-25/9))

    # horizontal jet
    j = plume(ex, SimpleJet(release_angle=0.0, k2=a, k3=b))
    @test isa(j, GasDispersion.SimpleJetSolution)
    @test isa(j, Plume)
    @test j(x,y,z) ≈ c
    @test j(-x,y,z) == 0.0

end
