@testset "Britter-McQuaid plume correlations" begin

    @test GasDispersion._bm_pl_c(-0.75)[2] ≈ [1.75, 1.92, 2.08, 2.25, 2.4, 2.6]
    @test GasDispersion._bm_pl_c(-0.40)[2] ≈ [1.7839999999999998, 2.016, 2.21, 2.3939999999999997, 2.564, 2.714]
    @test GasDispersion._bm_pl_c(-0.20)[2] ≈ [1.8319999999999999, 2.06, 2.25, 2.45, 2.63, 2.77]
    @test GasDispersion._bm_pl_c(0.0)[2] ≈ [1.78, 1.96, 2.16, 2.35, 2.56, 2.71]

    # test domain warning
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pl_c_1(1.1)
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pl_c_05(1.1)
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pl_c_02(1.1)
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pl_c_01(1.1)
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pl_c_005(1.1)
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pl_c_002(1.1)


end

@testset "Britter-McQuaid plume tests" begin
    # Britter-McQuaid example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 122
    R = 8.31446261815324
    ṁ = 0.23*425.6 # liquid spill rate times liquid density
    ρ = 1.76
    T = 273.15-162
    P = 101325
    MW = ρ*R*T/P
    Q = ṁ/ρ # gas volumetric flowrate: mass flowrate divided by gas density
    A = (π/4)*(1)^2 # release area, assuming a diameter of 1m.
    u = Q/A

    s = Substance(name = :LNG,
                  molar_weight=MW,
                  vapor_pressure=nothing,
                  gas_density = ρ,
                  liquid_density = 425.6,
                  reference_temp=T,
                  reference_pressure=P,
                  k = 0,
                  boiling_temp = 111.6, # Methane, NIST Webbook
                  latent_heat = 509880,  # J/kg, Methane
                  gas_heat_capacity = 2240, # J/kg/K, Methane
                  liquid_heat_capacity = 3349) # J/kg/K, Methane

    r = HorizontalJet( mass_rate = ṁ,
                  duration = ρ,
                  diameter = 1,
                  velocity = u,
                  height = 0,
                  pressure = P,
                  temperature = T,
                  fraction_liquid = 0)
    a = SimpleAtmosphere(windspeed=10.9, temperature=298, stability=ClassF)
    scn = Scenario(s,r,a)
    # known answers
    # initial concentration
    c₀ = 1
    # in the near-field
    x₁, c₁ = 2.26*20, 0.6720234544594696
    # example, in the interpolation region
    x₂, c₂ = 367, 0.050717667650511944
    # far field
    x₃, c₃ = 1200, 0.004921409666659286

    # test overall solution
    pl = plume(scn, BritterMcQuaidPlume)
    @test isa(pl, GasDispersion.BritterMcQuaidPlumeSolution)
    @test isa(pl, Plume)
    @test pl(x₁,0,0) ≈ c₁
    @test pl(x₂,0,0) ≈ c₂
    @test pl(x₃,0,0) ≈ c₃

    # test plume extent
    Lu  = 0.5*pl.D + 2*pl.lb
    @test pl(-Lu - 2*eps(Float64), 0, 0) == 0.0
    @test pl(-Lu + 2*eps(Float64), 0, 0) ≈ c₀

    Lh0 = pl.D + 8*pl.lb
    @test pl(0, Lh0 + 2*eps(Float64), 0) == 0.0
    @test pl(0, Lh0 - 2*eps(Float64), 0) ≈ c₀
    @test pl(0 - 2*eps(Float64), Lh0 + 2*eps(Float64), 0) == 0.0
    @test pl(0 - 2*eps(Float64), Lh0 - 2*eps(Float64), 0) ≈ c₀

    Lv = pl.D^2/(2*Lh0)
    @test pl(0, 0, -1) == 0.0
    @test pl(0, 0, Lv + 2*eps(Float64)) == 0.0
    @test pl(0, 0, Lv - 2*eps(Float64)) ≈ c₀

    # test large α near-field
    a = SimpleAtmosphere(windspeed=0.8322, temperature=298, stability=ClassF)
    scn = Scenario(s,r,a)
    pl_nf = plume(scn, BritterMcQuaidPlume)
    c_nf  = 0.5311215219329275
    @test pl_nf(pl_nf.xnf*pl_nf.D,0,0) ≈ c_nf

end
