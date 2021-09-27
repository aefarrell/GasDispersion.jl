@testset "Gaussian plume tests" begin

    @test_throws ErrorException plume(ambient, model="gaussian")

    # mising params for plume rise
    s = Scenario(test_scenario; release_temperature=missing)
    @test_throws ErrorException plume(s, model="gaussian", plumerise=true)

    s = Scenario(test_scenario; windspeed=3.5, jet_velocity=1.5)
    no_downwash = plume(s, model="gaussian")
    downwash = plume(s, model="gaussian", downwash=true)
    @test no_downwash(10, 0, 1) > downwash(10,0,1)

    @testset "Gaussian plume tests for class $class" for class in ["A", "B", "C", "D", "E", "F"]
        s1 = Scenario( test_scenario; pasquill_gifford=class )

        s2 = Scenario( test_scenario; release_temperature=s1.ambient_temperature,
                                      pasquill_gifford=class )
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
        s = Scenario( test_scenario; pasquill_gifford=class )

        h = s.release_height
        u = s.windspeed
        x = 10
        t = x/u

        test_puff = puff(s, model="gaussian")
        @test test_puff(x,0,h,t) > 0

    end
end
