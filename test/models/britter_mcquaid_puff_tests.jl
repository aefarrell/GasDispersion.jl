@testset "Britter-McQuaid puff tests" begin
    sub = Substance(name="test gas",
                    gas_density=1.2,
                    liquid_density=1000.,
                    reference_temp=298.,
                    reference_pressure=101325.0,
                    boiling_temp=100.,
                    latent_heat=1.,
                    gas_heat_capacity=1.,
                    liquid_heat_capacity=1.)
    rel = Release(mass_rate=1.0,
                  duration=Inf,
                  diameter=1.0,
                  velocity=1.0,
                  height=1.0,
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0.0)
    atm = DryAir(temperature=298.0,
                 pressure=101325.0,
                 windspeed=2.0,
                 stability=ClassF)
    scn = Scenario(sub,rel,atm)
    @test_throws ErrorException puff(scn, BritterMcQuaidPuff)

end
