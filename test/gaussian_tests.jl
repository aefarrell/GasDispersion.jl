@testset "Gaussian plume tests" begin

    @test_throws ErrorException plume(ambient, model="gaussian")

    @testset "Gaussian plume tests for class $class" for class in ["A", "B", "C", "D", "E", "F"]
        s1 = Scenario(
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

        s2 = Scenario(
            test_scenario.mass_emission_rate,
            test_scenario.release_duration,
            test_scenario.jet_diameter,
            test_scenario.jet_velocity,
            test_scenario.jet_density,
            test_scenario.release_pressure,
            test_scenario.ambient_temperature, # release at ambient temperature
            test_scenario.release_height,
            test_scenario.windspeed,
            test_scenario.ambient_density,
            test_scenario.ambient_pressure,
            test_scenario.ambient_temperature,
            class,   # pasquill stability class
        )
        h = s1.release_height

        no_plume_rise_1 = plume(s1, model="gaussian")
        no_plume_rise_2 = plume(s2, model="gaussian")

        @test no_plume_rise_1(10,0,h) ≈ no_plume_rise_2(10,0,h)

        plume_rise_1 = plume(s1, model="gaussian", plumerise=true, downwash=true)
        @test no_plume_rise_1(10,0,h) > plume_rise_1(10,0,h)

        plume_rise_2 = plume(s2, model="gaussian", plumerise=true, downwash=true)
        @test no_plume_rise_2(10,0,h) > plume_rise_2(10,0,h)
    end

end

@testset "Gaussian puff tests" begin
    @test_throws ErrorException puff(ambient, model="gaussian")

    @testset "Gaussian puff tests for class $class" for class in ["A", "B", "C", "D", "E", "F"]
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

        h = s.release_height
        u = s.windspeed
        x = 10
        t = x/u

        test_puff = puff(s, model="gaussian")
        @test test_puff(x,0,h,t) > 0

    end
end
