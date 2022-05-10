@testset "Britter-McQuaid plume tests" begin
    # Britter-McQuaid example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 122
    r = Release( mass_rate = (0.23*425.6),
                 duration = Inf,
                 diameter = 0,
                 velocity = 0,
                 height = 10.0,
                 pressure = 0,
                 temperature = (273.15-162),
                 density = 1.76)
    a = Ambient(windspeed=10.9, density=1.224, temperature=298, stability="F")
    ex = Scenario(r,a)
    # known answers
    # set 1 is covers the short-distance correlation
    # set 2 is the answer from the above example, which uses the interpolations
    x₁, c₁ = 2.26*20, 1.1827612798486666
    x₂, c₂ = 367.0, 0.0872919843565787

    # test type inheritance
    @test isa(plume(ex, BritterMcQuaidPlume()), Plume)

    @testset "Britter-McQuaid plume tests for class $class" for class in ["A", "B", "C", "D", "E", "F"]
        # because the windspeed is at 10m, the class should not impact the
        # calculations but this is a check that getting the corresponding
        # correlation doesn't throw any errors or do anything deeply strange
        a = Ambient(windspeed=10.9, density=1.224, temperature=298, stability=class)
        s = Scenario(r,a)
        britter_mcquaid = plume(s, BritterMcQuaidPlume())
        @test britter_mcquaid(x₁,0,0) ≈ c₁
        @test britter_mcquaid(x₂,0,0) ≈ c₂

    end
end

@testset "Britter-McQuaid puff tests" begin
    r = Release( mass_rate = (0.23*425.6),
                 duration = Inf,
                 diameter = 0,
                 velocity = 0,
                 height = 10.0,
                 pressure = 0,
                 temperature = (273.15-162),
                 density = 1.76)
    a = Ambient(windspeed=10.9, density=1.224, temperature=298, stability="F")
    ex = Scenario(r,a)

    @test_throws ErrorException puff(ex, BritterMcQuaidPuff())

end
