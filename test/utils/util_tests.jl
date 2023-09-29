# this is probably a bad way of doing things but I want to test some utils
# that are not otherwise exposed to the user
include("../../src/utils/utils.jl")

@testset "Property getters" begin
    sub = Substance("test",8.90,490.0,298.15,120935.0368,100.,123456., 78910.,4.19)
    rel = Release(1.0, 10.0, 0.25, 15.67, 2.0, 101325.0, 450., 0.67)
    atm = DryAir(100e3,273.15,287.05,2.0,5.0,ClassA)
    scn = Scenario(sub,rel,atm)

    @test _atmosphere_temperature(scn) == _temperature(atm) == 273.15
    @test _release_temperature(scn) == _temperature(rel) == 450.0
    @test _atmosphere_pressure(scn) == _pressure(atm) == 100e3
    @test _release_pressure(scn) == _pressure(rel) == 101325.0

    @test _mass_rate(scn) == _mass_rate(rel) == 1.0
    @test _duration(scn) == _duration(rel) == 10.0
    @test _release_mass(scn) == _mass(rel) == 10.0
    @test _release_diameter(scn) == _diameter(rel) == 0.25
    @test _release_area(scn) == _area(rel) ≈ (π/4)*0.25^2
    @test _release_velocity(scn) == _velocity(rel) == 15.67
    @test _release_flowrate(scn) == _flowrate(rel) ≈ (π/4)*(0.25^2)*15.67
    @test _release_height(scn) == _height(rel) == 2.0
    @test _release_liquid_fraction(scn) == _liquid_fraction(rel) == 0.67

    @test _windspeed(scn) == _windspeed(atm) == _velocity(atm) == 2.0
    @test _windspeed_height(scn) == _windspeed_height(atm) == _height(atm) == 5.0
    @test _stability(scn) == _stability(atm) == ClassA

end

@testset "Density functions" begin
    sub1 = Substance("test",8.90,490.0,298.15,120935.0368,100.,123456., 78910.,4.19)
    sub2 = Substance("test",(x,y)->y*x^2,(x,y)->y*x^3,2,3,100.,123456., 78910.,4.19)
    rel = Release(1.0, 10.0, 0.25, 15.67, 2.0, 101325.0, 450., 0.5)
    atm = DryAir(100e3,273.15,287.05,2.0,5.0,ClassA)
    scn1 = Scenario(sub1,rel,atm)
    scn2 = Scenario(sub2,rel,atm)

    @test _liquid_density(sub1) ≈ 490.0
    @test _liquid_density(sub2) ≈ 24
    @test _gas_density(sub1) ≈ 8.90
    @test _gas_density(sub2) ≈ 12

    @test _density(sub1, 0.5, 298.15, 120935.0368) ≈ 2/(1/8.9 + 1/490.)
    @test _density(sub2, 0.5, 2, 3) ≈ 48/3

    @test _release_density(scn1) == _density(sub1,0.5,450.0,101325.0)
    @test _release_density(scn2) == _density(sub2,0.5,450.0,101325.0)

    @test _density(atm, 100, 200) ≈ 2/287.05
    @test _atmosphere_density(scn1) == _density(atm)

end

@testset "Monin-Obukhov length tests" begin
    @test _monin_obukhov(1.2, ClassA) ≈ -11.609752888076077
    @test _monin_obukhov(1.2, ClassB) ≈ -26.818480014823884
    @test _monin_obukhov(1.2, ClassC) ≈ -129.91505611802876
    @test _monin_obukhov(1.2, ClassD) == Inf
    @test _monin_obukhov(1.2, ClassE) ≈ 129.91505611802876
    @test _monin_obukhov(1.2, ClassF) ≈ 26.818480014823884
end

@testset "Pasquill-Gifford dispersion tests" begin

    # Plume dispersion
    knowns = [(ClassA, 0.49842921341962687, 79.47756697338727),
              (ClassB, 0.3688140515374544, 0.15901388834172905),
              (ClassC, 0.247447127229602, 0.11437251727122057),
              (ClassD, 0.16025147287250416, 0.05115043235859766),
              (ClassE, 0.12018860465437811, 0.028796954615478678),
              (ClassF, 0.0794187446441675, 0.014462956106986533)]
    @testset "Stability class $class" for (class,cwind,vert) in knowns
        @test crosswind_dispersion(1.2, Plume, class, DefaultSet()) ≈ cwind
        @test vertical_dispersion(1.2, Plume, class, DefaultSet()) ≈ vert
    end

    # Puff dispersion
    knowns = [(ClassA, 0.21287234847910436, 0.6879188103871441),
              (ClassB, 0.16556738215041453, 0.6054498545741843),
              (ClassC, 0.11826241582172466, 0.3869880921561032),
              (ClassD, 0.07095744949303479, 0.17041904657983328),
              (ClassE, 0.047304966328689864, 0.11258170198247626),
              (ClassF, 0.023523465599668385, 0.05588182287353654)]
    @testset "Stability class $class" for (class,cwind,vert) in knowns
        @test crosswind_dispersion(1.2, Puff, class, DefaultSet()) ≈ cwind
        @test downwind_dispersion(1.2, Puff, class, DefaultSet()) ≈ cwind
        @test vertical_dispersion(1.2, Puff, class, DefaultSet()) ≈ vert
    end
end

@testset "Windspeed by powerlaw" begin

    u0, z0, p = 3.0, 1.0, 0.108
    a = DryAir(windspeed=u0, windspeed_height=z0, stability=ClassA)
    s = Scenario(Substance(:null,0,0,0,0,0,0,0,0),Release(0,0,0,0,1.0,0,0,0),a)
    @test _windspeed(s) == _windspeed(a) ≈ u0
    @test _windspeed(s,10) == _windspeed(a,10) == _windspeed(u0,z0,10,ClassA,DefaultSet())

    knowns = [(ClassA, 3.846991747968065),
              (ClassB, 3.882587524349958),
              (ClassC, 3.954770215669221),
              (ClassD, 4.160267486615665),
              (ClassE, 4.7876374417101975),
              (ClassF, 5.371817562105883)]

    @testset "Stability class $class" for (class, ans) in knowns
        @test  _windspeed(u0,z0,10,class,DefaultSet()) ≈ ans
    end

end

@testset "Windspeed by Monin-Obukhov length" begin
    ustar, zR, k = 3.0, 1.0, 0.35
    λ = _monin_obukhov(zR, ClassA)
    a(z) = (1-15*(z/λ))^0.25
    Ψ(a) = 2*log((1+a)/2) + log((1+a^2)/2) - 2*atan(a) + π/2

    @test _windspeed(10, ustar, zR, λ, ClassA) ≈ (ustar/k)*(log((10+zR)/zR) - Ψ(a(10)))
    @test _windspeed(10, ustar, zR, λ, ClassD) ≈ (ustar/k)*(log((10+zR)/zR))
    @test _windspeed(10, ustar, zR, λ, ClassE) ≈ (ustar/k)*(log((10+zR)/zR) - 4.7*(10/λ))
end

@testset "Plume rise" begin

    # test null case
    @test isa(NoPlumeRise(),PlumeRise)

    # test Briggs model of plume rise
    g = 9.80616 # m/s^2
    Tₐ = 288.15 # ambient temperature, K
    u = 4. # windspeed, m/s
    x = 50. # test point, m

    # Buoyant unstable plume
    # Fb = 50 case
    Tᵣ = 400 # K
    uⱼ = 72.93819699672669 # m/s
    xf = 565.0050541491846 # m
    Δhf = 100.71365158671999 # m
    sln = plume_rise(1.0,uⱼ,Tᵣ,u,Tₐ,ClassA)
    @test isa(sln,BuoyantPlume)
    @test sln ≈ BuoyantPlume(50.0,xf,u,Δhf)
    @test plume_rise(x, sln) ≈ 20.0
    @test plume_rise(2xf, sln) ≈ Δhf

    # Fb = 60 case
    uⱼ = 87.52583639607202
    xf = 612.07897481385
    Δhf = 112.8895989623143
    sln = plume_rise(1.0,uⱼ,Tᵣ,u,Tₐ,ClassA)
    @test isa(sln,BuoyantPlume)
    @test sln ≈ BuoyantPlume(60.0,xf,u,Δhf)

    # Momentum dominated unstable plume
    # Fb = 50 case
    Tᵣ = 325 # K
    uⱼ = 179.87751923862007 # m/s
    xf = 565.0050541491846 # m
    Fm = 7171.814541070672
    β = (1/3) + (u/uⱼ)
    Δhf = 3*(uⱼ/u) # m
    sln = plume_rise(1.0,uⱼ,Tᵣ,u,Tₐ,ClassA)
    @test isa(sln,MomentumPlume)
    @test sln ≈ MomentumPlume(Fm,xf,β,nothing,u,Δhf,ClassA)
    @test plume_rise(x, sln) ≈ 81.01824072514351
    @test plume_rise(2xf, sln) ≈ Δhf

    # Buoyant stable plume
    # Fb = 50, class E
    Tᵣ = 400 # K
    uⱼ = 72.93819699672669 # m/s
    xf = 317.60677324769046 # m
    Δhf = 68.59722859012221 # m
    sln = plume_rise(1.0,uⱼ,Tᵣ,u,Tₐ,ClassE)
    @test isa(sln,BuoyantPlume)
    @test sln ≈ BuoyantPlume(50.0,xf,u,Δhf)

    # Fb = 50, class F
    Tᵣ = 400 # K
    uⱼ = 72.93819699672669 # m/s
    xf = 240.0881533494489 # m
    Δhf = 56.92380039947288 # m
    sln = plume_rise(1.0,uⱼ,Tᵣ,u,Tₐ,ClassF)
    @test isa(sln,BuoyantPlume)
    @test sln ≈ BuoyantPlume(50.0,xf,u,Δhf)

    # Momentum dominated stable plume
    # Fb = 67, class E
    Tᵣ = 325 # K
    uⱼ = 240.91531595745576 # m/s
    Fm = 12864.831225945449
    xf = 240.83782417699823 # m
    β = (1/3)+(u/uⱼ)
    s = 0.0006806288391462781
    Δhf = 19.04522445600697 # m
    sln = plume_rise(1.0,uⱼ,Tᵣ,u,Tₐ,ClassE)
    @test isa(sln,MomentumPlume)
    @test sln ≈ MomentumPlume(Fm,xf,β,s,u,Δhf,ClassE)
    @test plume_rise(x, sln) ≈ Δhf
    @test plume_rise(2xf, sln) ≈ Δhf

    # Momentum dominated stable plume
    # Fb = 67, class F
    Tᵣ = 325 # K
    uⱼ = 240.91531595745576 # m/s
    Fm = 12864.831225945449
    xf = 182.0562825914961 # m
    β = (1/3)+(u/uⱼ)
    s = 0.0011911004685059867
    Δhf = 17.349211998937378 # m
    sln = plume_rise(1.0,uⱼ,Tᵣ,u,Tₐ,ClassF)
    @test isa(sln,MomentumPlume)
    @test sln ≈ MomentumPlume(Fm,xf,β,s,u,Δhf,ClassF)

end

@testset "Equation Sets" begin
    
    include("ccps_tests.jl")

end