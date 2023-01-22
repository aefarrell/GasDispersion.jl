@testset "Britter-McQuaid plume tests" begin
    # Britter-McQuaid example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 122
    s = Substance(name = :test,
                  gas_density = 1.76,
                  liquid_density = 425.6,
                  reference_temp=(273.15-162),
                  reference_pressure=101325.0,
                  boiling_temp = 111.6, # Methane, NIST Webbook
                  latent_heat = 8.17/0.0160425, # Methane, NIST Webbook
                  gas_heat_capacity = 1.6151, # Methane, NIST Webbook
                  liquid_heat_capacity = 2.0564) # Methane, NIST Webbook
    r = Release( mass_rate = (0.23*425.6),
                 duration = 174,
                 diameter = 0,
                 velocity = 0,
                 height = 10.0,
                 pressure = 101325.0,
                 temperature = (273.15-162),
                 fraction_liquid = 0.0)
    a = DryAir(windspeed=10.9, temperature=298, stability=ClassF)
    scn = Scenario(s,r,a)
    # known answers
    # set 1 is covers the short-distance correlation
    # set 2 is the answer from the above example, which uses the interpolations
    x₁, c₁ = 2.26*20, 1.1827612798486664
    x₂, c₂ = 367.0, 0.08922951841378725

    # test type inheritance
    pl = plume(scn, BritterMcQuaidPlume)
    @test isa(pl, GasDispersion.BritterMcQuaidPlumeSolution)
    @test isa(pl, Plume)
    @test pl(x₁,0,0) ≈ c₁
    @test pl(x₂,0,0) ≈ c₂

    # test plume extent
    Lu  = 0.5*pl.D + 2*pl.lb
    @test pl(-Lu - 2*eps(Float64), 0, 0) == 0.0

    Lh0 = pl.D + 8*pl.lb
    @test pl(0, Lh0 + 2*eps(Float64), 0) == 0.0
    @test pl(0 - 2*eps(Float64), Lh0 + 2*eps(Float64), 0) == 0.0

    Lv = pl.D^2/Lh0
    @test pl(0, 0, -1) == 0.0
    @test pl(0, 0, Lv + 2*eps(Float64)) == 0.0

end
