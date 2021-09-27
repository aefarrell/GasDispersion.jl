@testset "Britter-McQuaid plume tests" begin

    @test_throws ErrorException plume(ambient, model="britter-mcquaid")

    @testset "Britter-McQuaid plume tests for class $class" for class in ["A", "B", "C", "D", "E", "F"]
        s = Scenario(
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
            class,   # pasquill stability class
        )
        britter_mcquaid = plume(s, model="britter-mcquaid")
        @test britter_mcquaid(10,0,0) > 0.0

    end
end
