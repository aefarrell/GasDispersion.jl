@testset "Britter-McQuaid puff correlations" begin

    @test GasDispersion._bm_pf_c(-0.75)[2] ≈ [0.7, 0.85, 0.95, 1.15, 1.48, 1.83, 2.075]
    @test GasDispersion._bm_pf_c(0.0)[2] ≈ [0.81, 1.0, 1.19, 1.39, 1.62, 1.83, 2.05]
    @test GasDispersion._bm_pf_c(0.75)[2] ≈ [0.93, 1.03, 1.1849999999999998, 1.375, 1.525, 1.68, 1.8474999999999997]

    # test domain warning
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pf_c_1(1.1)
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pf_c_05(1.1)
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pf_c_02(1.1)
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pf_c_01(1.1)
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pf_c_005(1.1)
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pf_c_002(1.1)
    @test_logs (:warn, "α= 1.1 is out of range") GasDispersion._bm_pf_c_001(1.1)

end

@testset "Britter-McQuaid puff tests" begin
    # Britter-McQuaid example, *Workbook on the Dispersion of Dense Gases*
    # Potchefstroom Accident
    s = Substance(name = :pseudo_ammonia_air_mix,
                  molar_weight=0.017031,
                  vapor_pressure=0,
                  gas_density = 1.434,
                  liquid_density = 681.63,        # Ammonia, NIST Webbook
                  reference_temp=(273.15-33.316), # boiling point of Ammonia, NIST Webbook
                  reference_pressure=101325,
                  k=0,
                  boiling_temp = (273.15-33.316), # Ammonia, NIST Webbook
                  latent_heat = 8.17/0.0160425,   # Ammonia, NIST Webbook
                  gas_heat_capacity = 1.6151,     # Ammonia, NIST Webbook
                  liquid_heat_capacity = 2.0564)  # Ammonia, NIST Webbook
    r = HorizontalJet( mass_rate = 7.25e4*1.434,
                 duration = 1,
                 diameter = 1,
                 velocity = 7.25e4/(0.25*π),
                 height = 0,
                 pressure = 101325,
                 temperature = (273.15-33.316),   # Ammonia, NIST Webbook
                 fraction_liquid = 0)
    a = SimpleAtmosphere(windspeed=2, temperature=293.15, stability=ClassF())
    scn = Scenario(s,r,a)
    # known answers
    # initial concentration
    c₀ = 1/16.4
    # in the near-field
    x₁, t₁, c₁ = 160, 160/(0.4*2), 0.0030830245167018195
    # examples, in the interpolation region
    x₂, t₂, c₂ = 333, 124, 0.006097560975609757
    x₃, t₃, c₃ = 333, 1396, 0.00018497619507887978
    # example in the far-field
    x₄, t₄, c₄ = 2500, 2500/(0.4*2), 0.0007259996788424369/16.4

    # test overall solution
    pf = puff(scn, BritterMcQuaidPuff(); temp_correction=false)
    @test isa(pf, GasDispersion.BritterMcQuaidPuffSolution)
    @test isa(pf, Puff)
    @test c₀*pf(x₁,0,0,t₁) ≈ c₁
    @test c₀*pf(x₂,0,0,t₂) ≈ c₂
    @test c₀*pf(x₃,0,0,t₃) ≈ c₃
    @test c₀*pf(x₄,0,0,t₄) ≈ c₄

    # test puff extent
    R² = cbrt(3*pf.V₀/4π)^2 + 1.2*√(pf.gₒ*pf.V₀)*t₁
    @test c₀*pf(x₁, √(R²), 0, t₁) ≈ c₁
    @test pf(x₁, √(R²) + 1, 0, t₁) == 0.0

    H = (c₀*pf.V₀)/(c₁*π*R²)
    @test c₀*pf(x₁, 0, H, t₁) ≈ c₁
    @test pf(x₁, 0, H + 1, t₁) == 0.0

    # test temperature correction
    pf = puff(scn, BritterMcQuaidPuff(); temp_correction=true)
    @test pf.T′ ≈ (273.15-33.316) / 293.15

end
