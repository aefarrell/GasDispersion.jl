@testset "Britter-McQuaid puff tests" begin
    s = Substance(name = :test,
                  gas_density = 1.76,
                  liquid_density = 1000.0,
                  boiling_temp = NaN,
                  latent_heat = NaN,
                  gas_heat_capacity = NaN,
                  liquid_heat_capacity = NaN)
    r = Release( mass_rate = (0.23*425.6),
                 duration = Inf,
                 diameter = 0,
                 velocity = 0,
                 height = 10.0,
                 pressure = 0,
                 temperature = (273.15-162),
                 fraction_liquid = 0.0)
    a = Ambient(windspeed=10.9, density=1.224, temperature=298, stability=ClassF)
    ex = Scenario(s,r,a)

    @test_throws ErrorException puff(ex, BritterMcQuaidPuff)

end
