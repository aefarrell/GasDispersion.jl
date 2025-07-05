@testset "Gaussian mixing layer tests" begin
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
    hj = HorizontalJet(mass_rate=0.1,
                  duration=Inf,
                  diameter=1,
                  velocity=1,
                  height=0,
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0)

    hj_toohigh = HorizontalJet(mass_rate=0.1,
                  duration=Inf,
                  diameter=1,
                  velocity=1,
                  height=10_000,  # too high for the mixing layer
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0)
    
    vj = VerticalJet(mass_rate=0.1,
                  duration=Inf,
                  diameter=1,
                  velocity=1,
                  height=10,
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0)

    vj_toohigh = VerticalJet(mass_rate=0.1,
                  duration=Inf,
                  diameter=1,
                  velocity=1,
                  height=10_000,  # too high for the mixing layer
                  temperature=298.,
                  pressure=101325.,
                  fraction_liquid=0)
    
    vj_rh_high = VerticalJet(mass_rate=0.1,
                  duration=Inf,
                  diameter=1,
                  velocity=100,
                  height=630,  # plume rise will be too high
                  temperature=500.,
                  pressure=101325.,
                  fraction_liquid=0)

    stbl = SimpleAtmosphere(temperature=298,
                 pressure=101325,
                 windspeed=2,
                 windspeed_height=10,
                 stability=ClassF())

    neut = SimpleAtmosphere(temperature=298,
                 pressure=101325,
                 windspeed=2,
                 windspeed_height=10,
                 stability=ClassD())

    
    # horizontal jet test default behaviour and type inheritance
    @test plume(Scenario(sub,hj,stbl), GaussianMixingLayer()).verticalterm isa GasDispersion.SimpleVerticalTerm
    @test plume(Scenario(sub,hj,neut), GaussianMixingLayer()).verticalterm isa GasDispersion.SimpleMixingLayer
    @test plume(Scenario(sub,hj,neut), GaussianMixingLayer(); method=:periodicmixinglayer).verticalterm isa GasDispersion.PeriodicMixingLayer
    @test_throws ErrorException plume(Scenario(sub,hj,neut), GaussianMixingLayer(); method=:someothermethod)
    @test_throws ErrorException plume(Scenario(sub,hj_toohigh,neut), GaussianMixingLayer())
    
    # vertical jet test default behaviour and type inheritance
    @test plume(Scenario(sub,vj,stbl), GaussianMixingLayer()).verticalterm isa GasDispersion.SimpleVerticalTerm
    @test plume(Scenario(sub,vj,neut), GaussianMixingLayer()).verticalterm isa GasDispersion.SimpleMixingLayer
    @test plume(Scenario(sub,vj,neut), GaussianMixingLayer(); method=:periodicmixinglayer).verticalterm isa GasDispersion.PeriodicMixingLayer
    @test plume(Scenario(sub,vj,neut), GaussianMixingLayer(); plumerise=false).plumerise isa GasDispersion.NoPlumeRise
    @test_throws ErrorException plume(Scenario(sub,vj,neut), GaussianMixingLayer(); method=:someothermethod)
    @test_throws ErrorException plume(Scenario(sub,vj_toohigh,neut), GaussianMixingLayer())

    # testing limit on plume rise
    pl_rise = plume(Scenario(sub,vj_rh_high,neut), GaussianMixingLayer(); downwash=false, plumerise=true)
    @test pl_rise.plumerise.final_rise ≈ pl_rise.verticalterm.mixing_height
    pl_low_terms = plume(Scenario(sub,vj_rh_high,neut), GaussianMixingLayer(); downwash=false, plumerise=true, n_terms=1)
    @test pl_low_terms(pl_rise.plumerise.xf,0,0.9*pl_rise.verticalterm.mixing_height) < pl_rise(pl_rise.plumerise.xf,0,0.9*pl_rise.verticalterm.mixing_height)
    
    # integration test with a scenario -- simple mixing layer
    pl_smpl = plume(Scenario(sub,vj,neut), GaussianMixingLayer(); downwash=true, plumerise=true)
    @test pl_smpl.verticalterm.mixing_height ≈ 640.0311954829524
    @test pl_smpl(500, 0, 0) ≈ 1.7194314353343084e-5

    # integration test with a scenario -- periodic mixing layer
    pl_per = plume(Scenario(sub,vj,neut), GaussianMixingLayer(); downwash=true, plumerise=true, method=:periodicmixinglayer, n_terms=100_000_000)
    @test pl_per.verticalterm.mixing_height ≈ 640.0311954829524
    @test pl_per(500, 0, 0) ≈ 1.7194314353343084e-5

end
