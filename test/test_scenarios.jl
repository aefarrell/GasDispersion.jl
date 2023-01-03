# test constructor
test_scenario = Scenario(Substance("test",8.90,490.0,298.15,120935.0368,0,0,0,0),
                 Release(1.0, 10.0, 0.25, 15.67, 1.0, 101325.0, 450, 1.0),
                 Ambient())
@test isa(test_scenario, Scenario)
@test test_scenario ≈ test_scenario

# I have shamelessly stolen this from the tests for show()
replstr(x, kv::Pair...) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80), kv...), MIME("text/plain"), x), x)
test_scenario_str = "Substance: test \nRelease conditions:\n    mass_rate: 1.0 kg/s \n    duration: 10.0 s \n    diameter: 0.25 m \n    velocity: 15.67 m/s \n    height: 1.0 m \n    pressure: 101325.0 Pa \n    temperature: 450 K \n    fraction_liquid: 1.0  \nAtmospheric conditions:\n    pressure: 101325 Pa \n    temperature: 298.15 K \n    density: 1.225 kg/m^3 \n    windspeed: 1.5 m/s \n    windspeed_height: 10 m \n    stability: ClassF  \n"
@test replstr(test_scenario) == test_scenario_str

@testset "Scenario Builder tests" begin
    # Liquid jet example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 40
    a = Ambient()
    s1 = Substance(name="test liquid", gas_density=0.0, liquid_density=490.0,
         boiling_temp=0.0,latent_heat=0.0,gas_heat_capacity=0.0,
         liquid_heat_capacity=0.0)
    #known answer
    ljet = Scenario(s1, Release(mass_rate=0.21691154763598,duration=Inf,
            diameter=0.01,velocity=5.636333880812954,height=1.0,pressure=101325,
            temperature=298.15,fraction_liquid=1.0), a)
    @test ljet ≈ scenario_builder(s1,JetSource,a;phase=:liquid,dischargecoef=0.63,diameter=0.01,pressure=120935.0368,temperature=298.15,height=1.0)

    # Gas jet example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 47
    s2 = Substance(name="test gas",gas_density=8.90,liquid_density=0.0,
         reference_temp=298.0,reference_pressure=501e3,boiling_temp=0.0,
         latent_heat=0.0,gas_heat_capacity=0.0,liquid_heat_capacity=0.0)
    #known answer
    gjet = Scenario(s2, Release(mass_rate=0.09002799947040846,duration=Inf,
            diameter=0.01,velocity=208.58711308961637,height=1.0,
            pressure=287766.01316878956,temperature=277.2093023255814,
            fraction_liquid = 0.0),a)
    @test gjet ≈ scenario_builder(s2,JetSource,a;phase=:gas,dischargecoef=0.85,k=1.15,diameter=0.01,pressure=501e3,temperature=298.0,height=1.0)

    # testing invalid phase
    @test_throws ErrorException scenario_builder(s1,JetSource;phase=:fake,
                        diameter=0,pressure=0,temperature=0,height=0)

end
