@testset "Gaussian plume tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 97
    sub = Substance(name="test gas",
                    gas_density=1.2268,
                    liquid_density=1000.,
                    reference_temp=298.,
                    reference_pressure=101325.,
                    boiling_temp=100.,
                    latent_heat=1.,
                    gas_heat_capacity=1.,
                    liquid_heat_capacity=1.)
    rel = Release(mass_rate=0.1,
                  duration=Inf,
                  diameter=1.0,
                  velocity=1.0,
                  height=0.0,
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0.0)
    atm = DryAir(temperature=298.0,
                 pressure=101325.0,
                 windspeed=2.0,
                 stability=ClassF)
    scn = Scenario(sub,rel,atm)
    pl = plume(scn)

    # test default behaviour and type inheritance
    @test isa(pl, GasDispersion.GaussianPlumeSolution)
    @test isa(pl, Plume)
    @test pl(-1,0,0) == 0.0
    @test pl(0,0,0) ≈ 0.1/(π/4)
    @test pl(500,0,0) ≈ 0.00010346728324507407

    # test stack downwash calculation
    rel = Release(mass_rate=0.1,
                  duration=Inf,
                  diameter=1.0,
                  velocity=1.0,
                  height=10.0,
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0.0)
    scn = Scenario(sub,rel,atm)
    pl = plume(scn; downwash=true)
    @test pl.effective_stack_height ≈ 8.0

    # test with plume rise
    pl = plume(scn; plumerise=true)
    @test isa(pl.plumerise,GasDispersion.BriggsModel)
    @test pl(500, 0, 10) ≈ 5.297761895423843e-5

end
