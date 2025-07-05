@testset "Depreciated behaviour tests" begin
    # an example scenario that should work with all models
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
                  duration=3600,
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
                 stability=ClassF())
    scn = Scenario(sub,rel,atm)

    @test_logs (:warn, "Instantiating SimpleAtmosphere with a StabilityClass type is deprecated, use an instance instead. For example, use stability=ClassF() instead of stability=ClassF.") begin
        SimpleAtmosphere(temperature=298, pressure=101325, windspeed=2, windspeed_height=1, stability=ClassF)
    end

    @test_logs (:warn, "scenario_builder(substance, JetSource, atmosphere) is deprecated, use scenario_builder(substance, JetSource(), atmosphere) instead.") begin
        scenario_builder(sub, JetSource, atm; phase=:gas,dischargecoef=0.85,diameter=0.01,pressure=501e3,temperature=298,height=1)
    end

    @test_logs (:warn, "plume(scenario, GaussianPlume, eqs) is deprecated, use plume(scenario, GaussianPlume(), eqs) instead.") begin
        plume(scn, GaussianPlume, DefaultSet)
    end

    @test_logs (:warn, "plume(scenario, SimpleJet, eqs) is deprecated, use plume(scenario, SimpleJet(), eqs) instead.") begin
        plume(scn, SimpleJet, DefaultSet)
    end

    @test_logs (:warn, "plume(scenario, BritterMcQuaidPlume, eqs) is deprecated, use plume(scenario, BritterMcQuaidPlume(), eqs) instead.") begin
        plume(scn, BritterMcQuaidPlume, DefaultSet)
    end

    @test_logs (:warn, "puff(scenario, GaussianPuff, eqs) is deprecated, use puff(scenario, GaussianPuff(), eqs) instead.") begin
        puff(scn, GaussianPuff, DefaultSet)
    end

    @test_logs (:warn, "puff(scenario, IntPuff, eqs) is deprecated, use puff(scenario, IntPuff(), eqs) instead.") begin
        puff(scn, IntPuff, DefaultPuffSet)
    end

    @test_logs (:warn, "puff(scenario, BritterMcQuaidPuff, eqs) is deprecated, use puff(scenario, BritterMcQuaidPuff(), eqs) instead.") begin
        puff(scn, BritterMcQuaidPuff, DefaultSet)
    end

    @test_logs (:warn, "puff(scenario, SLAB, eqs) is deprecated, use puff(scenario, SLAB(), eqs) instead.") begin
        puff(scn, SLAB, DefaultSet)
    end


end