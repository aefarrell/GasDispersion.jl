include("../src/utils/pasquill_gifford.jl")
include("../src/models/plume_rise.jl")

@testset "Gaussian plume tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 97
    ex = Scenario( Dict([
        :mass_emission_rate => 0.1,
        :jet_diameter => 10.0,
        :jet_velocity => 1.0,
        :release_height => 0.0,
        :windspeed => 2.0,
        :pasquill_gifford => "F"
    ]))
    # known answers
    x₁ = 500.0

    # test type inheritance
    @test isa(plume(ex, GaussianPlume()), Plume)

    # missing model params
    @test_throws MissingException plume(ambient, GaussianPlume())

    # missing params exclusively for plume rise
    s = Scenario(ex; release_temperature=missing)
    @test_throws MissingException plume(s, GaussianPlume(plumerise=true))

    # invalid stability class for plume rise
    s = Scenario(ex; release_temperature=325, ambient_temperature=298,
                     pasquill_gifford="not a real one")
    @test_throws ErrorException plume_rise(s, true)

    @testset "Gaussian plume tests for class $class" for class in ["A", "B", "C", "D", "E", "F"]
        s = Scenario(ex; pasquill_gifford=class)

        m, u = s.mass_emission_rate, s.windspeed
        σy, σz = crosswind_dispersion(class)(x₁), vertical_dispersion(class)(x₁)
        c₁ = m / (π*σy*σz*u)

        # basic model, no plume rise no stack downwash
        no_plume_rise = plume(s, GaussianPlume())
        @test no_plume_rise(x₁,0,0) ≈ c₁

        # basic model, with stack downwash included
        Dⱼ, uⱼ = s.jet_diameter, s.jet_velocity
        u, h = 2*uⱼ, 1.0
        hₑ = h + 2*Dⱼ*( (uⱼ/u) - 1.5 )
        c₂ = m / (2π*σy*σz*u)*( exp(-0.5*((h-hₑ)/σz)^2) + exp(-0.5*((h+hₑ)/σz)^2) )
        downwash = plume(Scenario(s; release_height=h, windspeed=u), GaussianPlume(downwash=true))
        @test downwash(x₁, 0, h) ≈ c₂

        # models with plume rise
        # buoyancy driven, Fb<55
        u = s.windspeed
        Tᵣ, Tₐ = 315, 298
        bs1 = Scenario(s; release_height=h, release_temperature=Tᵣ,
                                            ambient_temperature=Tₐ)
        Δh = plume_rise(bs1, true)(x₁)
        hₑ  = h + Δh
        σyₑ = √( (Δh/3.5)^2 + σy^2 )
        σzₑ = √( (Δh/3.5)^2 + σz^2 )
        c₃ = m / (2π*σyₑ*σzₑ*u)*( exp(-0.5*((h-hₑ)/σzₑ)^2) + exp(-0.5*((h+hₑ)/σzₑ)^2) )
        b1 = plume(Scenario(bs1), GaussianPlume(plumerise=true))
        @test b1(x₁, 0, h) ≈ c₃

        # buoyancy driven, Fb>55
        Tᵣ = 400
        bs2 = Scenario(bs1; release_temperature=Tᵣ)
        Δh = plume_rise(bs2, true)(x₁)
        hₑ  = h + Δh
        σyₑ = √( (Δh/3.5)^2 + σy^2 )
        σzₑ = √( (Δh/3.5)^2 + σz^2 )
        c₄ = m / (2π*σyₑ*σzₑ*u)*( exp(-0.5*((h-hₑ)/σzₑ)^2) + exp(-0.5*((h+hₑ)/σzₑ)^2) )
        b2 = plume(Scenario(bs2), GaussianPlume(plumerise=true))
        @test b2(x₁, 0, h) ≈ c₄

        # momentum driven, Fb<55
        Tᵣ = Tₐ
        ms1 = Scenario(bs1; release_temperature=Tᵣ)
        Δh = plume_rise(ms1, true)(x₁)
        hₑ  = h + Δh
        σyₑ = √( (Δh/3.5)^2 + σy^2 )
        σzₑ = √( (Δh/3.5)^2 + σz^2 )
        c₅ = m / (2π*σyₑ*σzₑ*u)*( exp(-0.5*((h-hₑ)/σzₑ)^2) + exp(-0.5*((h+hₑ)/σzₑ)^2) )
        m1 = plume(Scenario(ms1), GaussianPlume(plumerise=true))
        @test m1(x₁, 0, h) ≈ c₅

        # momentum driven, Fb>55
        g = 9.80616
        uⱼ = ( (4/0.00575)*55*g )^(3/5) * (1/Dⱼ)
        ms2 = Scenario(bs1; jet_velocity=uⱼ)
        Δh = plume_rise(ms2, true)(x₁)
        hₑ  = h + Δh
        σyₑ = √( (Δh/3.5)^2 + σy^2 )
        σzₑ = √( (Δh/3.5)^2 + σz^2 )
        c₆ = m / (2π*σyₑ*σzₑ*u)*( exp(-0.5*((h-hₑ)/σzₑ)^2) + exp(-0.5*((h+hₑ)/σzₑ)^2) )
        m2 = plume(Scenario(ms2), GaussianPlume(plumerise=true))
        @test m2(x₁, 0, h) ≈ c₆

    end

end

@testset "Gaussian puff tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 97
    ex = Scenario( Dict([
        :mass_emission_rate => 1.0,
        :release_duration => 1.0,
        :jet_diameter => 10.0,
        :jet_velocity => 1.0,
        :release_height => 0.0,
        :windspeed => 2.0,
        :pasquill_gifford => "F"
    ]))
    # knowns
    m = ex.mass_emission_rate*ex.release_duration
    h = ex.release_height
    u = ex.windspeed
    x₁ = 500.0
    t₁ = x₁/u

    # test type inheritance
    @test isa(puff(ex, GaussianPuff()), Puff)

    # missing model params
    @test_throws MissingException puff(ambient, GaussianPuff())

    @testset "Gaussian puff tests for class $class" for class in ["A", "B", "C", "D", "E", "F"]
        s = Scenario( ex; pasquill_gifford=class )

        σx = crosswind_dispersion(class, plume=false)(x₁)
        σy = σx
        σz = vertical_dispersion(class, plume=false)(x₁)

        c₁ = m / (√(2π)*π*σx*σy*σz)

        test_puff = puff(s, GaussianPuff())
        @test test_puff(x₁,0,h,t₁) ≈ c₁

    end
end
