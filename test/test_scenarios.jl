s = Substance("test",8.90,490.0,0,0,0,0)
r = Release(1.0, 10.0, 0.25, 15.67, 1.0, 101325.0, 450, 1.0)
a = Ambient()
test_scenario = Scenario(s,r,a)

# test constructor
@test isa(test_scenario, Scenario)
@test test_scenario ≈ Scenario(s,r,a)

# I have shamelessly stolen this from the tests for show()
replstr(x, kv::Pair...) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80), kv...), MIME("text/plain"), x), x)
test_scenario_str = "Atmospheric conditions:\n    pressure: 101325 Pa \n    temperature: 298.15 K \n    density: 1.225 kg/m^3 \n    windspeed: 1.5 m/s \n    stability: F  \n    \nRelease conditions:\n    mass_rate: 1.0 kg/s \n    duration: 10.0 s \n    diameter: 0.25 m \n    velocity: 15.67 m/s \n    height: 1.0 m \n    pressure: 101325.0 Pa \n    temperature: 450 K \n    fraction_liquid: 1.0  \n    "

@test replstr(test_scenario) == test_scenario_str

@testset "Scenario Builder tests" begin
    # Liquid jet example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 40
    js1 = JetSource(phase=:liquid,
            dischargecoef=0.63,
            diameter=0.01,
            pressure=120935.0368,
            temperature=298.15,
            height=1.0 )
    a = Ambient()
    ljet = Scenario(s,Release( mass_rate = 0.21691154763598,
                duration = Inf,
                diameter = 0.01,
                velocity = 5.636333880812954,
                height = 1.0,
                pressure = 101325,
                temperature = 298.15,
                fraction_liquid = 1.0), a )
    @test ljet == scenario_builder(s,js1,a)

    # Gas jet example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 47
    js2 = JetSource( phase=:gas,
            dischargecoef=0.85,
            k=1.15,
            diameter=0.01,
            pressure=501e3,
            temperature=298.15,
            height=1.0)
    gjet = Scenario(Substance{Float64, Int64}(
                        "test", 5.495411838308172, 490.0, 0, 0, 0, 0),
                Release( mass_rate = 0.09002799947040846,
                        duration = Inf,
                        diameter = 0.01,
                        velocity = 208.58711308961637,
                        height = 1.0,
                        pressure = 287766.01316878956,
                        temperature = 277.3488372093023,
                        fraction_liquid = 0.0), a )
    @test gjet == scenario_builder(s,js2,a)

    # testing default behaviour
    @test scenario_builder(s,js1) == scenario_builder(s,js1, a)

    # testing invalid phase
    @test_throws ErrorException scenario_builder(s,JetSource( phase=:fake,
                        diameter=0,pressure=0,temperature=0,height=0))

end
