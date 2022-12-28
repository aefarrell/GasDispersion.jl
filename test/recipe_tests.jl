using Plots

@testset "Plot recipe tests" begin
    null_substance = Substance(name=:test,gas_density=NaN,liquid_density=NaN,
                    boiling_temp=NaN,latent_heat=NaN,gas_heat_capacity=NaN,
                    liquid_heat_capacity=NaN)
    ex = Scenario(null_substance,
                  Release(mass_rate=0.1,duration=10.0,diameter=10.0,
                   velocity=1.0,height=0.0,pressure=0.0,temperature=0.0,
                   fraction_liquid=0.0),
                  Ambient(windspeed=2.0, stability=:F))

    pl = plume(ex)
    p1 = plot(pl; height=0, xlims=(-1,100), ylims=(-50,50), clims=(0,pl.max_concentration))
    @test isa(p1,Plots.Plot)

    pf = puff(ex)
    p2 = plot(pf, 10; height=0, xlims=(17.5,22.5), ylims=(-2.5,2.5))
    @test isa(p2,Plots.Plot)

end
