@testset "Gaussian plume tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 97
    sub = Substance(name="test gas",
                    molar_weight=1,
                    vapor_pressure=nothing,
                    gas_density=1.2268,
                    liquid_density=1000,
                    reference_temp=298,
                    reference_pressure=101325,
                    k=0,
                    boiling_temp=100,
                    latent_heat=1,
                    gas_heat_capacity=1,
                    liquid_heat_capacity=1)
    rel = HorizontalJet(mass_rate=0.1,
                  duration=Inf,
                  diameter=1,
                  velocity=1,
                  height=0,
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0)
    atm = SimpleAtmosphere(temperature=298,
                 pressure=101325,
                 windspeed=2,
                 windspeed_height=1,
                 stability=ClassF)
    scn = Scenario(sub,rel,atm)
    pl = plume(scn)

    # test default behaviour and type inheritance
    @test GasDispersion.GaussianPlumeSolution(scn,:test,1.0,2,3,4,5,GasDispersion.SimpleCrossTerm(),GasDispersion.SimpleVerticalTerm(),GasDispersion.NoPlumeRise(),ClassA,DefaultSet(),GasDispersion.ProblemDomain(0.0,Inf,-Inf,Inf,0.0,Inf)) isa GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.SimpleCrossTerm, GasDispersion.SimpleVerticalTerm, GasDispersion.NoPlumeRise, ClassA, GasDispersion.BasicEquationSet{GasDispersion.DefaultWind, Nothing, GasDispersion.Defaultσy, GasDispersion.Defaultσz}, GasDispersion.ProblemDomain{Float64}}
    @test pl isa GasDispersion.GaussianPlumeSolution
    @test pl isa Plume
    @test pl(-1,0,0) == 0.0
    @test pl(0,0,0) ≈ 0.1/(π/4)/1.2268
    @test pl(500,0,0) ≈ 0.00010346728324507407/1.2268

    # test stack downwash calculation
    rel = VerticalJet(mass_rate=0.1,
                  duration=Inf,
                  diameter=1,
                  velocity=1,
                  height=10,
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0)
    atm = SimpleAtmosphere(temperature=298,
                  pressure=101325,
                  windspeed=2,
                  windspeed_height=10,
                  stability=ClassF)
    scn = Scenario(sub,rel,atm)
    pl = plume(scn; downwash=true, plumerise=false)
    @test pl.effective_stack_height ≈ 8

    # test with plume rise
    pl = plume(scn; plumerise=true)
    @test pl.plumerise isa GasDispersion.BriggsModel
    @test pl(-1,0,0) == 0.0
    @test pl(0,0,10) ≈ 0.1/(π/4)/1.2268
    @test pl(500, 0, 10) ≈ 5.297761895423843e-5/1.2268

end
