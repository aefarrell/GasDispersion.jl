include("../../src/utils/pasquill_gifford.jl")
include("../../src/utils/plume_rise.jl")

null_substance = Substance(name=:test,gas_density=NaN,liquid_density=NaN,
                boiling_temp=NaN,latent_heat=NaN,gas_heat_capacity=NaN,
                liquid_heat_capacity=NaN)

@testset "Gaussian plume tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 97

    ex = Scenario(null_substance,
                  Release(mass_rate=0.1,duration=Inf,diameter=10.0,
                    velocity=1.0,height=0.0,pressure=0.0,temperature=0.0,
                    fraction_liquid=0.0),
                  Ambient(windspeed=2.0, stability=ClassF))
    x₁ = 500.0
    pl = plume(ex)

    # test default behaviour and type inheritance
    @test isa(pl,GasDispersion.GaussianPlumeSolution)
    @test isa(pl, Plume)
    @test pl(-1,0,0) == 0.0
    @test pl(0,0,0) == pl.max_concentration

    @testset "Gaussian plume tests for class $class" for class in [ClassA, ClassB, ClassC, ClassD, ClassE, ClassF]

        s = Scenario(null_substance,
                      Release(mass_rate=0.1,duration=Inf,diameter=10.0,
                        velocity=1.0,height=0.0,pressure=0.0,temperature=0.0,
                        fraction_liquid=0.0),
                      Ambient(windspeed=2.0, stability=class))

        m, u = _mass_rate(s), _windspeed(s)
        σy = crosswind_dispersion(x₁, Plume, class)
        σz = vertical_dispersion(x₁, Plume, class)
        c₁ = m / (π*σy*σz*u)

        # basic model, no plume rise no stack downwash
        no_plume_rise = plume(s, GaussianPlume)
        @test no_plume_rise(x₁,0,0) ≈ c₁

        # basic model, with stack downwash included
        Dⱼ, uⱼ = _release_diameter(s), _release_velocity(s)
        u, h = 2*uⱼ, 1.0
        hₑ = h + 2*Dⱼ*( (uⱼ/u) - 1.5 )
        c₂ = m / (2π*σy*σz*u)*( exp(-0.5*((h-hₑ)/σz)^2) + exp(-0.5*((h+hₑ)/σz)^2) )

        dw = Scenario(null_substance,
                    Release(mass_rate=m, duration=Inf, diameter=Dⱼ,
                     velocity=uⱼ, height=h, pressure=0, temperature=0,
                     fraction_liquid=0),
                    Ambient(windspeed=u, stability=class))
        downwash = plume(dw, GaussianPlume; downwash=true)
        @test downwash(x₁, 0, h) ≈ c₂
    end

end

@testset "Gaussian puff tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 97
    ex = Scenario(null_substance,
                Release(mass_rate = 0.1, duration = 10.0, diameter = 10.0,
                velocity = 1.0, height = 0.0, pressure = 0, temperature = 0,
                fraction_liquid = 0),
                Ambient(windspeed=2.0, stability=ClassF))
    # knowns
    m = _release_mass(ex)
    h = _release_height(ex)
    u = _windspeed(ex)
    x₁ = 500.0
    t₁ = x₁/u
    pf = puff(ex)

    # test default behaviour and type inheritance
    @test isa(pf,GasDispersion.GaussianPuffSolution)
    @test isa(pf, Puff)
    @test pf(-1,0,h,t₁) == 0.0

    @testset "Gaussian puff tests for class $class" for class in [ClassA, ClassB, ClassC, ClassD, ClassE, ClassF]
        s = Scenario(null_substance,
                    Release(mass_rate = 0.1, duration = 10.0, diameter = 10.0,
                    velocity = 1.0, height = 0.0, pressure = 0, temperature = 0,
                    fraction_liquid = 0),
                    Ambient(windspeed=2.0, stability=class))

        σx = crosswind_dispersion(x₁, Puff, class)
        σy = σx
        σz = vertical_dispersion(x₁, Puff, class)

        c₁ = m / (√(2π)*π*σx*σy*σz)

        test_puff = puff(s, GaussianPuff)
        @test test_puff(x₁,0,h,t₁) ≈ c₁

    end
end

@testset "Integrated Gaussian puff tests" begin
    # Gaussian plume example, *Guidelines for Consequence Analysis of Chemical
    # Releases* CCPS, 1999, pg 97
    ex = Scenario(null_substance,
                Release(mass_rate = 0.1, duration = 10.0, diameter = 10.0,
                velocity = 1.0, height = 0.0, pressure = 0, temperature = 0,
                fraction_liquid = 0),
                Ambient(windspeed=2.0, stability=ClassF))
    h = _release_height(ex)
    u = _windspeed(ex)
    x₁ = 500.0
    t₁ = x₁/u
    Δt = _duration(ex)

    # testing default behaviour
    @test isa(puff(ex, IntPuff;n=1), GasDispersion.GaussianPuffSolution)
    @test isa(puff(ex, IntPuff;n=3), GasDispersion.IntPuffSolution{<:Integer,<:StabilityClass})
    @test isa(puff(ex, IntPuff), GasDispersion.IntPuffSolution{<:Float64,<:StabilityClass})
    @test_throws ErrorException puff(ex, IntPuff; n=0)

    # testing 3 puffs
    gp = puff(ex, GaussianPuff)
    ip = puff(ex, IntPuff;n=3)
    @test ip(x₁,0,h,t₁) ≈ (1/3)*(gp(x₁,0,h,t₁) + gp(x₁,0,h,t₁-0.5*Δt) + gp(x₁,0,h,t₁-Δt))
    @test ip(-1,0,h,t₁) == 0.0

    # testing ∞ puffs
    ip∞ = puff(ex,IntPuff)
    @test ip∞(x₁,0,h,t₁) ≈ 0.000711709646491573
    @test ip∞(-1,0,h,t₁) == 0.0
end
