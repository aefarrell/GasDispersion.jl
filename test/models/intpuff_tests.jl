@testset "Integrated Gaussian puff tests" begin
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
                 stability=ClassF())
    scn = Scenario(sub,rel,atm)
    x₁, t₁, Δt, h = 500, 250, 10, 10

    # testing default behaviour
    @test GasDispersion.IntPuffSolution(scn,:test_promotion,1.0,2,3,4,5,6,DefaultPuffSet) isa GasDispersion.IntPuffSolution{Float64, Int64, GasDispersion.BasicEquationSet{GasDispersion.DefaultWind, GasDispersion.CCPSPuffσx, GasDispersion.CCPSPuffσy, GasDispersion.CCPSPuffσz}}
    @test puff(scn, IntPuff();n=1) isa GasDispersion.GaussianPuffSolution
    @test puff(scn, IntPuff();n=3) isa GasDispersion.IntPuffSolution{<:Number,<:Integer,<:EquationSet}
    @test puff(scn, IntPuff()) isa GasDispersion.PalazziSolution{<:Number,<:EquationSet}
    @test_throws ErrorException puff(scn, IntPuff(); n=0)

    # testing 3 puffs
    gp = puff(scn, GaussianPuff())
    ip = puff(scn, IntPuff();n=3)
    @test ip(x₁,0,h,t₁) ≈ (1/3)*(gp(x₁,0,h,t₁) + gp(x₁,0,h,t₁-0.5*Δt) + gp(x₁,0,h,t₁-Δt))
    @test ip(x₁,0,h,-t₁) == 0.0

    # testing ∞ puffs
    ip∞ = puff(scn,IntPuff())
    @test ip∞(x₁,0,h,t₁) ≈ 0.00035586710616021256/1.2268
    @test ip∞(x₁,0,h,-t₁) == 0.0
end
