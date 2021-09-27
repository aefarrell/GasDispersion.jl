ambient = Scenario(
    missing,   # mass emission rate, kg/s
    missing,  # release duration, s
    missing,  # jet diameter, m
    missing, # jet velocity, m/s
    missing,   # jet density, kg/m^3
    101325.0,# release_pressure, Pa
    298.15,   # release temperature, K
    missing,   # release height, m
    missing,   # windspeed, m/s
    1.225, # ambient density, kg/m^3
    101325.0,# ambient pressure, Pa
    298.15,# ambient temperature, K
    missing,   # pasquill stability class
)

test_scenario = Scenario(
    1.0,   # mass emission rate, kg/s
    10.0,  # release duration, s
    0.25,  # jet diameter, m
    15.67, # jet velocity, m/s
    1.3,   # jet density, kg/m^3
    101325.0,# release_pressure, Pa
    450,   # release temperature, K
    1.0,   # release height, m
    1.5,   # windspeed, m/s
    1.225, # ambient density, kg/m^3
    101325.0,# ambient pressure, Pa
    298.15,# ambient temperature, K
    "F",   # pasquill stability class
)

bad_class = Scenario(
    test_scenario.mass_emission_rate,
    test_scenario.release_duration,
    test_scenario.jet_diameter,
    test_scenario.jet_velocity,
    test_scenario.jet_density,
    test_scenario.release_pressure,
    test_scenario.release_temperature,
    test_scenario.release_height,
    test_scenario.windspeed,
    test_scenario.ambient_density,
    test_scenario.ambient_pressure,
    test_scenario.ambient_temperature,
    "error",
)



@testset "Scenario constructor tests" begin
    d = Dict{Symbol,Union{Missing, Number, String}}([
        :release_pressure => 101325.0,   # Pa
        :release_temperature => 298.15,# K
        :ambient_density => 1.225,     # kg/m^3
        :ambient_pressure => 101325.0,   # Pa
        :ambient_temperature => 298.15,# K
    ])

    @test Scenario(d) == ambient

    @test Scenario(test_scenario) == test_scenario
    @test Scenario(test_scenario; pasquill_gifford="error") == bad_class

end

@testset "Scenario Builder tests" begin

    jet = Dict([
        :mass_emission_rate => 0.21691154763598,
        :jet_diameter => 0.01,
        :jet_velocity => 5.636333880812954,
        :jet_density => 490.0,
        :release_pressure => 120935.0368,
        :release_temperature => 298.15,
        :release_height => 0.0,
        :windspeed => 1.5,
        :ambient_density => 1.225,
        :ambient_pressure => 101325.0,
        :ambient_temperature => 298.15,
        :pasquill_gifford => "F"
    ])

    @test Scenario(jet) == scenario_builder(jet[:release_pressure], jet[:release_temperature];
                           model="jet", phase="liquid",
                           stability=jet[:pasquill_gifford], windspeed=jet[:windspeed],
                           liquid_density=jet[:jet_density], hole_diameter=jet[:jet_diameter],
                           discharge_coeff=0.63)

    @test_throws ErrorException scenario_builder(jet[:release_pressure], jet[:release_temperature];
                           model="something else", phase="liquid",
                           stability=jet[:pasquill_gifford], windspeed=jet[:windspeed])

end
