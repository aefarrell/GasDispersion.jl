test_scenario = Scenario(Release(1.0, 10.0, 0.25, 15.67, 1.0, 101325.0, 450, 1.3),
                         Ambient() )

# I have shamelessly stolen this from the tests for show()
replstr(x, kv::Pair...) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80), kv...), MIME("text/plain"), x), x)
test_scenario_str = "Atmospheric conditions:\n    pressure: 101325 Pa \n    temperature: 298.15 K \n    density: 1.225 kg/m^3 \n    windspeed: 1.5 m/s \n    stability: F  \n    \nRelease conditions:\n    mass_rate: 1.0 kg/s \n    duration: 10.0 s \n    diameter: 0.25 m \n    velocity: 15.67 m/s \n    height: 1.0 m \n    pressure: 101325.0 Pa \n    temperature: 450 K \n    density: 1.3 kg/m^3 \n    "

@test replstr(test_scenario) == test_scenario_str

@testset "Scenario Builder tests" begin
    # Liquid jet example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 40
    js = JetSource( phase=:liquid,
            dischargecoef=0.63,
            diameter=0.01,
            pressure=120935.0368,
            temperature=298.15,
            density=490.0,
            height=1.0 )
    a = Ambient()
    jet = Scenario( Release( mass_rate = 0.21691154763598,
                duration = Inf,
                diameter = 0.01,
                velocity = 5.636333880812954,
                height = 1.0,
                pressure = 120935.0368,
                temperature = 298.15,
                density = 490.0 ), a )
    # using default ambient properties
    @test jet == scenario_builder(js,a)


end
