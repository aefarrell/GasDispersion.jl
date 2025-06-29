@testset "Gaussian mixing layer tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 97
    sub = Substance(name="test gas",
                    molar_weight=1,
                    vapor_pressure=nothing,
                    gas_density=1.2268,
                    liquid_density=1000,
                    reference_temp=298,
                    reference_pressure=101325,
                    k=0,
                    boiling_temp=100,
                    latent_heat=1,
                    gas_heat_capacity=1,
                    liquid_heat_capacity=1)
    rel = HorizontalJet(mass_rate=0.1,
                  duration=Inf,
                  diameter=1,
                  velocity=1,
                  height=0,
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0)
    atm1 = SimpleAtmosphere(temperature=298,
                 pressure=101325,
                 windspeed=2,
                 windspeed_height=1,
                 stability=ClassF)

    atm2 = SimpleAtmosphere(temperature=298,
                 pressure=101325,
                 windspeed=2,
                 windspeed_height=1,
                 stability=ClassD)
    scn = Scenario(sub,rel,atm1)
    
    # test default behaviour and type inheritance
    pl1 = plume(scn, GaussianMixingLayer)
    @test pl1.verticalterm isa GasDispersion.SimpleVerticalTerm
    
    scn2 = Scenario(sub,rel,atm2)
    pl2 = plume(scn2, GaussianMixingLayer)
    @test pl2.verticalterm isa GasDispersion.SimpleMixingLayer

    # test vertical jet with mixing layer
    rel = VerticalJet(mass_rate=0.1,
                  duration=Inf,
                  diameter=1,
                  velocity=1,
                  height=10,
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0)
    scn3 = Scenario(sub,rel,atm1)
    pl3 = plume(scn3, GaussianMixingLayer)
    @test pl3.verticalterm isa GasDispersion.SimpleVerticalTerm
    
    atm3 = SimpleAtmosphere(temperature=298,
                  pressure=101325,
                  windspeed=2,
                  windspeed_height=10,
                  stability=ClassD)
    scn4 = Scenario(sub,rel,atm3)
    pl4 = plume(scn4, GaussianMixingLayer; downwash=true, plumerise=true)
    @test pl4.verticalterm isa GasDispersion.SimpleMixingLayer
    @test pl4.verticalterm.mixing_height ≈ 640.0311954829524
    @test pl4(500, 0, 0) ≈ 1.7194314353343084e-5

end
