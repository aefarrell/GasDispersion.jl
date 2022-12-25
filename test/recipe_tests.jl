using Plots

@testset "Plot recipe tests" begin
    ex = Scenario(Release(mass_rate = 0.1, duration = Inf, diameter = 10.0,
                velocity = 1.0, height = 0.0, pressure = 0, temperature = 0,
                density = 0), Ambient(windspeed=2.0, stability="F"))

    p = plume(ex)
    @test isa(plot(p),Plots.Plot)

end
