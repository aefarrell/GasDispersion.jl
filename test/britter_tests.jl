@testset "Britter-McQuaid plume tests" begin
    # Britter-McQuaid example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 122
    ex = Scenario( Dict([
    :mass_emission_rate => (0.23*425.6),
    :release_height => 10.0,
    :jet_density => 1.76,
    :release_temperature => (273.15-162),
    :windspeed => 10.9,
    :ambient_density => 1.224,
    :ambient_temperature => 298,
    :pasquill_gifford => "F"
    ]))
    # known answers
    # set 1 is covers the short-distance correlation
    # set 2 is the answer from the above example, which uses the interpolations
    x₁, c₁ = 2.26*20, 1.1827612798486666
    x₂, c₂ = 367.0, 0.0872919843565787

    # missing model parameters
    @test_throws MissingException plume(ambient, model=:brittermcquaid)

    # test type inheritance
    @test isa(plume(ex, model=:brittermcquaid), PlumeModel)

    @testset "Britter-McQuaid plume tests for class $class" for class in ["A", "B", "C", "D", "E", "F"]
        # because the windspeed is at 10m, the class should not impact the
        # calculations but this is a check that getting the corresponding
        # correlation doesn't throw any errors or do anything deeply strange
        s = Scenario( ex; pasquill_gifford=class )
        britter_mcquaid = plume(s, model=:brittermcquaid)
        @test britter_mcquaid(x₁,0,0) ≈ c₁
        @test britter_mcquaid(x₂,0,0) ≈ c₂

    end
end

@testset "Britter-McQuaid puff tests" begin
    ex = Scenario( Dict([
    :mass_emission_rate => (0.23*425.6),
    :release_height => 10.0,
    :jet_density => 1.76,
    :release_temperature => (273.15-162),
    :windspeed => 10.9,
    :ambient_density => 1.224,
    :ambient_temperature => 298,
    :pasquill_gifford => "F"
    ]))

    @test_throws ErrorException puff(ex, model=:brittermcquaid)

end
