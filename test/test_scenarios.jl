test_scenario = Scenario(
    1.0,   # mass emission rate, kg/s
    10.0,  # release duration, s
    0.25,  # jet diameter, m
    15.67, # jet velocity, m/s
    1.3,   # jet density, kg/m^3
    101325,# release_pressure, Pa
    450,   # release temperature, K
    1.0,   # release height, m
    1.5,   # windspeed, m/s
    1.225, # ambient density, kg/m^3
    101325,# ambient pressure, Pa
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
