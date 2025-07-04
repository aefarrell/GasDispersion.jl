# the unit tests are inspired by those used by Measurements.jl

using RecipesBase

# this is a hack, but it doesn't work otherwise ¯\_(ツ)_/¯
RecipesBase.is_key_supported(k::Symbol) = true

@testset "Plot recipe tests" begin
    # default scenario to test
    sub = Substance(name="test gas",
                    molar_weight=1.0,
                    vapor_pressure=1.0,
                    gas_density=1.2268,
                    liquid_density=1000.,
                    reference_temp=298.,
                    reference_pressure=101325.,
                    boiling_temp=100.,
                    latent_heat=1.,
                    gas_heat_capacity=1.,
                    liquid_heat_capacity=1.)
    rel = HorizontalJet(mass_rate=0.1,
                  duration=10,
                  diameter=1.0,
                  velocity=1.0,
                  height=0.0,
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0.0)
    atm = SimpleAtmosphere(temperature=298.0,
                 pressure=101325.0,
                 windspeed=2.0,
                 stability=ClassF())
    scn = Scenario(sub,rel,atm)

    # expected recipe data
    res = Dict{Symbol,Any}(
            :ylims       => (-50, 50),
            :xlims       => (-1, 1000),
            :ylabel      => "Crosswind distance, m",
            :ny          => 100,
            :nx          => 100,
            :height      => 2,
            :seriestype  => :contour,
            :seriescolor => :thermal,
            :fill        => true,
            :xlabel      => "Downwind distance, m")

    # testing plume
    pl = plume(scn)
    rec, = RecipesBase.apply_recipe(Dict{Symbol,Any}(),pl)
    @test getfield(rec,1) == res
    @test rec.args[1] == range(-1,1000; length=100)
    @test rec.args[2] == range(-50,50; length=100)
    @test rec.args[3](10,10) == pl(10,10,2)

    # testing puff
    pf = puff(scn)
    rec, = RecipesBase.apply_recipe(Dict{Symbol,Any}(),pf)
    @test getfield(rec,1) == res
    @test rec.args[1] == range(-1,1000; length=100)
    @test rec.args[2] == range(-50,50; length=100)
    @test rec.args[3](10,10) == pf(10,10,2,0)

end
