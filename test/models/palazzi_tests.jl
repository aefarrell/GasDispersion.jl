@testset "Palazzi short duration puff tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 101
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
                  duration=10,
                  diameter=1,
                  velocity=1,
                  height=10,
                  temperature=298,
                  pressure=101325,
                  fraction_liquid=0)
    atm = SimpleAtmosphere(temperature=298,
                 pressure=101325,
                 windspeed=2,
                 stability=ClassF)
    scn = Scenario(sub,rel,atm)
    x₁, t₁, Δt, h = 500, 250, 10, 10

    # testing default behaviour
    @test puff(scn, Palazzi) isa GasDispersion.PalazziSolution{<:Number,<:StabilityClass,<:EquationSet}
    pf = puff(scn, Palazzi)
    @test pf(x₁,0,h,t₁) ≈ 0.00035586710616021256/1.2268
    @test pf(x₁,0,h,-t₁) == 0.0
    
end
