using Plots

@testset "Plot recipe tests" begin
    ex = Scenario(Release(mass_rate = 0.1, duration = 10.0, diameter = 10.0,
                velocity = 1.0, height = 0.0, pressure = 0, temperature = 0,
                density = 0), Ambient(windspeed=2.0, stability="F"))

    pl = plume(ex)
    p1 = plot(pl; height=0, xlims=(-1,100), ylims=(-50,50), clims=(0,pl.max_concentration))
    @test isa(p1,Plots.Plot)

    pf = puff(ex)
    p2 = plot(pf, 10; height=0, xlims=(17.5,22.5), ylims=(-2.5,2.5))
    @test isa(p2,Plots.Plot)

end
