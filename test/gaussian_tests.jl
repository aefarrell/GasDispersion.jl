include("../src/utils/pasquill_gifford.jl")
include("../src/models/plume_rise.jl")

@testset "Gaussian plume tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 97
    ex = Scenario(Release(mass_rate = 0.1, duration = Inf, diameter = 10.0,
                velocity = 1.0, height = 0.0, pressure = 0, temperature = 0,
                density = 0), Ambient(windspeed=2.0, stability="F"))
    # known answers
    x₁ = 500.0

    # test default behaviour
    @test plume(ex) == plume(ex, GaussianPlume())

    # test type inheritance
    @test isa(plume(ex, GaussianPlume()), Plume)

    # test for bad classes
    bad = Scenario(Release(mass_rate = 0.1, duration = Inf, diameter = 10.0,
                velocity = 1.0, height = 0.0, pressure = 0, temperature = 0,
                density = 0), Ambient(windspeed=2.0, stability="error"))
    @test_throws ErrorException plume(bad, GaussianPlume(plumerise=false))
    @test_throws ErrorException plume(bad, GaussianPlume(plumerise=true))

    @testset "Gaussian plume tests for class $class" for class in ["A", "B", "C", "D", "E", "F"]

        s = Scenario( Release(mass_rate = 0.1, duration = Inf, diameter = 10.0,
                    velocity = 1.0, height = 0.0, pressure = 0, temperature = 0,
                    density = 0), Ambient(windspeed=2.0, stability=class))

        m, u = s.release.mass_rate, s.atmosphere.windspeed
        σy, σz = crosswind_dispersion(class)(x₁), vertical_dispersion(class)(x₁)
        c₁ = m / (π*σy*σz*u)

        # basic model, no plume rise no stack downwash
        no_plume_rise = plume(s, GaussianPlume())
        @test no_plume_rise(x₁,0,0) ≈ c₁

        # basic model, with stack downwash included
        Dⱼ, uⱼ = s.release.diameter, s.release.velocity
        u, h = 2*uⱼ, 1.0
        hₑ = h + 2*Dⱼ*( (uⱼ/u) - 1.5 )
        c₂ = m / (2π*σy*σz*u)*( exp(-0.5*((h-hₑ)/σz)^2) + exp(-0.5*((h+hₑ)/σz)^2) )

        dw = Scenario( Release(mass_rate=m, duration=Inf, diameter=Dⱼ,
                     velocity=uⱼ, height=h, pressure=0, temperature=0,
                     density=0), Ambient(windspeed=u, stability=class))
        downwash = plume(dw, GaussianPlume(downwash=true))
        @test downwash(x₁, 0, h) ≈ c₂

        # models with plume rise
        # buoyancy driven, Fb<55
        u = s.atmosphere.windspeed
        Tᵣ, Tₐ = 315, 298
        bs1 = Scenario(Release(mass_rate=m, duration=Inf, diameter=Dⱼ,
                     velocity=uⱼ, height=h, pressure=0, temperature=Tᵣ,
                     density=0), Ambient(windspeed=u, stability=class, temperature=Tₐ))
        Δh = plume_rise(bs1, true)(x₁)
        hₑ  = h + Δh
        σyₑ = √( (Δh/3.5)^2 + σy^2 )
        σzₑ = √( (Δh/3.5)^2 + σz^2 )
        c₃ = m / (2π*σyₑ*σzₑ*u)*( exp(-0.5*((h-hₑ)/σzₑ)^2) + exp(-0.5*((h+hₑ)/σzₑ)^2) )
        b1 = plume(bs1, GaussianPlume(plumerise=true))
        @test b1(x₁, 0, h) ≈ c₃

        # buoyancy driven, Fb>55
        Tᵣ = 400
        bs2 = Scenario(Release(mass_rate=m, duration=Inf, diameter=Dⱼ,
                     velocity=uⱼ, height=h, pressure=0, temperature=Tᵣ,
                     density=0), Ambient(windspeed=u, stability=class, temperature=Tₐ))
        Δh = plume_rise(bs2, true)(x₁)
        hₑ  = h + Δh
        σyₑ = √( (Δh/3.5)^2 + σy^2 )
        σzₑ = √( (Δh/3.5)^2 + σz^2 )
        c₄ = m / (2π*σyₑ*σzₑ*u)*( exp(-0.5*((h-hₑ)/σzₑ)^2) + exp(-0.5*((h+hₑ)/σzₑ)^2) )
        b2 = plume(bs2, GaussianPlume(plumerise=true))
        @test b2(x₁, 0, h) ≈ c₄

        # momentum driven, Fb<55
        Tᵣ = Tₐ
        ms1 = Scenario(Release(mass_rate=m, duration=Inf, diameter=Dⱼ,
                     velocity=uⱼ, height=h, pressure=0, temperature=Tᵣ,
                     density=0), Ambient(windspeed=u, stability=class, temperature=Tₐ))
        Δh = plume_rise(ms1, true)(x₁)
        hₑ  = h + Δh
        σyₑ = √( (Δh/3.5)^2 + σy^2 )
        σzₑ = √( (Δh/3.5)^2 + σz^2 )
        c₅ = m / (2π*σyₑ*σzₑ*u)*( exp(-0.5*((h-hₑ)/σzₑ)^2) + exp(-0.5*((h+hₑ)/σzₑ)^2) )
        m1 = plume(ms1, GaussianPlume(plumerise=true))
        @test m1(x₁, 0, h) ≈ c₅

        # momentum driven, Fb>55
        Tᵣ = 315
        g  = 9.80616
        uⱼ = ( (4/0.00575)*55*g )^(3/5) * (1/Dⱼ)
        ms2 = Scenario(Release(mass_rate=m, duration=Inf, diameter=Dⱼ,
                     velocity=uⱼ, height=h, pressure=0, temperature=Tᵣ,
                     density=0), Ambient(windspeed=u, stability=class, temperature=Tₐ))
        Δh = plume_rise(ms2, true)(x₁)
        hₑ  = h + Δh
        σyₑ = √( (Δh/3.5)^2 + σy^2 )
        σzₑ = √( (Δh/3.5)^2 + σz^2 )
        c₆ = m / (2π*σyₑ*σzₑ*u)*( exp(-0.5*((h-hₑ)/σzₑ)^2) + exp(-0.5*((h+hₑ)/σzₑ)^2) )
        m2 = plume(ms2, GaussianPlume(plumerise=true))
        @test m2(x₁, 0, h) ≈ c₆

    end

end

@testset "Gaussian puff tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 97
    ex = Scenario(Release(mass_rate = 0.1, duration = 10.0, diameter = 10.0,
                velocity = 1.0, height = 0.0, pressure = 0, temperature = 0,
                density = 0), Ambient(windspeed=2.0, stability="F"))
    # knowns
    m = ex.release.mass_rate*ex.release.duration
    h = ex.release.height
    u = ex.atmosphere.windspeed
    x₁ = 500.0
    t₁ = x₁/u

    # test default behaviour
    @test puff(ex) == puff(ex, GaussianPuff())

    # test type inheritance
    @test isa(puff(ex, GaussianPuff()), Puff)

    # test for bad classes
    bad = Scenario(Release(mass_rate = 0.1, duration = Inf, diameter = 10.0,
                velocity = 1.0, height = 0.0, pressure = 0, temperature = 0,
                density = 0), Ambient(windspeed=2.0, stability="error"))
    @test_throws ErrorException puff(bad, GaussianPuff())

    @testset "Gaussian puff tests for class $class" for class in ["A", "B", "C", "D", "E", "F"]
        s = Scenario(Release(mass_rate = 0.1, duration = 10.0, diameter = 10.0,
                    velocity = 1.0, height = 0.0, pressure = 0, temperature = 0,
                    density = 0), Ambient(windspeed=2.0, stability=class))

        σx = crosswind_dispersion(class, plume=false)(x₁)
        σy = σx
        σz = vertical_dispersion(class, plume=false)(x₁)

        c₁ = m / (√(2π)*π*σx*σy*σz)

        test_puff = puff(s, GaussianPuff())
        @test test_puff(x₁,0,h,t₁) ≈ c₁

    end
end
