@testset "Gaussian puff tests" begin
    # Gaussian puff example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 101
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
                  duration=10,
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
    pf = puff(scn)
    # knowns
    x₁ = 500.0
    t₁ = 250.0
    σy = 5.047929716825321 # There is a mistake in the ref, it gives this as 6.1
    σz = 2.214836678536473
    c₁ = 1/(√(2π)*π*σy^2*σz)

    # test default behaviour and type inheritance
    @test isa(pf,GasDispersion.GaussianPuffSolution)
    @test isa(pf, Puff)
    @test pf(x₁,0,0,-t₁) == 0.0
    @test pf(x₁,0,0,t₁) ≈ c₁

end