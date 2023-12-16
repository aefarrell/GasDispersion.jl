@testset "Scenario Builder tests" begin
    # Liquid jet example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 40
    a = SimpleAtmosphere()
    s1 = Substance(name="test liquid", gas_density=0.0, liquid_density=490.0,
         boiling_temp=0.0,latent_heat=0.0,gas_heat_capacity=0.0,
         liquid_heat_capacity=0.0)
    #known answer
    ljet = Scenario(s1, HorizontalJet(mass_rate=0.21691154763598,duration=Inf,
            diameter=0.01,velocity=5.636333880812954,height=1.0,pressure=101325,
            temperature=298.15,fraction_liquid=1.0), a)
    @test ljet ≈ scenario_builder(s1,JetSource,a;phase=:liquid,dischargecoef=0.63,diameter=0.01,pressure=120935.0368,temperature=298.15,height=1.0)

    # Gas jet example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 47
    s2 = Substance(name="test gas",gas_density=8.90,liquid_density=0.0,
         reference_temp=298.0,reference_pressure=501e3,boiling_temp=0.0,
         latent_heat=0.0,gas_heat_capacity=0.0,liquid_heat_capacity=0.0)
    #known answer
    gjet = Scenario(s2, HorizontalJet(mass_rate=0.09002799947040846,duration=Inf,
            diameter=0.01,velocity=208.58711308961637,height=1.0,
            pressure=287766.01316878956,temperature=277.2093023255814,
            fraction_liquid = 0.0),a)
    @test gjet ≈ scenario_builder(s2,JetSource,a;phase=:gas,dischargecoef=0.85,k=1.15,diameter=0.01,pressure=501e3,temperature=298.0,height=1.0)

    # testing invalid phase
    @test_throws ErrorException scenario_builder(s1,JetSource;phase=:fake,
                        diameter=0,pressure=0,temperature=0,height=0)

    # testing default behaviour
    @test ljet ≈ scenario_builder(s1,JetSource;phase=:liquid,dischargecoef=0.63,diameter=0.01,pressure=120935.0368,temperature=298.15,height=1.0)
    
end
