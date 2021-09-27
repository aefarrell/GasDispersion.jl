@testset "Britter-McQuaid plume tests" begin

    @test_throws ErrorException plume(ambient, model="britter-mcquaid")

    @testset "Britter-McQuaid plume tests for class $class" for class in ["A", "B", "C", "D", "E", "F"]
        s = Scenario( test_scenario; pasquill_gifford=class )
        britter_mcquaid = plume(s, model="britter-mcquaid")
        # I should really add real values here
        x′_crit = 30*test_scenario.jet_diameter
        @test britter_mcquaid(0.5*x′_crit,0,0) > 0.0
        @test britter_mcquaid(2*x′_crit,0,0) > 0.0

    end
end
