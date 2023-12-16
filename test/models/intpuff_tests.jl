@testset "Integrated Gaussian puff tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
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
    rel = HorizontalJet(mass_rate=0.1,
                  duration=10.0,
                  diameter=1.0,
                  velocity=1.0,
                  height=10.0,
                  temperature=298.0,
                  pressure=101325.0,
                  fraction_liquid=0.0)
    atm = SimpleAtmosphere(temperature=298.0,
                 pressure=101325.0,
                 windspeed=2.0,
                 stability=ClassF)
    scn = Scenario(sub,rel,atm)
    x₁, t₁, Δt, h = 500.0, 250.0, 10.0, 10.0

    # testing default behaviour
    @test isa(puff(scn, IntPuff;n=1), GasDispersion.GaussianPuffSolution)
    @test isa(puff(scn, IntPuff;n=3), GasDispersion.IntPuffSolution{<:Integer,<:StabilityClass})
    @test isa(puff(scn, IntPuff), GasDispersion.IntPuffSolution{<:Float64,<:StabilityClass})
    @test_throws ErrorException puff(scn, IntPuff; n=0)

    # testing 3 puffs
    gp = puff(scn, GaussianPuff)
    ip = puff(scn, IntPuff;n=3)
    @test ip(x₁,0,h,t₁) ≈ (1/3)*(gp(x₁,0,h,t₁) + gp(x₁,0,h,t₁-0.5*Δt) + gp(x₁,0,h,t₁-Δt))
    @test ip(x₁,0,h,-t₁) == 0.0

    # testing ∞ puffs
    ip∞ = puff(scn,IntPuff)
    @test ip∞(x₁,0,h,t₁) ≈ 0.00035586710616021256/1.2268
    @test ip∞(x₁,0,h,-t₁) == 0.0
end
