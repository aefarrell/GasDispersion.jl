r = Release(
        1.0,   # mass emission rate, kg/s
        10.0,  # release duration, s
        0.25,  # jet diameter, m
        15.67, # jet velocity, m/s
        1.0,   # release height, m
        101325.0,# release_pressure, Pa
        450,   # release temperature, K
        1.3,   # jet density, kg/m^3
    )

test_scenario = Scenario(r, Ambient() )

bad_class = Scenario(r, Ambient(stability="error") )

# I have shamelessly stolen this from the tests for show()
replstr(x, kv::Pair...) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80), kv...), MIME("text/plain"), x), x)
test_scenario_str = "Atmospheric conditions:\n    pressure: 101325 Pa \n    temperature: 298.15 K \n    density: 1.225 kg/m^3 \n    windspeed: 1.5 m/s \n    stability: F  \n    \nRelease conditions:\n    mass_rate: 1.0 kg/s \n    duration: 10.0 s \n    diameter: 0.25 m \n    velocity: 15.67 m/s \n    height: 1.0 m \n    pressure: 101325.0 Pa \n    temperature: 450 K \n    density: 1.3 kg/m^3 \n    "
@testset "Scenario constructor tests" begin
    @test replstr(test_scenario) == test_scenario_str
end

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
