
@testset "Monin-Obukhov length tests" begin
    @test GasDispersion._monin_obukhov(1.2, ClassA) ≈ -11.609752888076077
    @test GasDispersion._monin_obukhov(1.2, ClassB) ≈ -26.818480014823884
    @test GasDispersion._monin_obukhov(1.2, ClassC) ≈ -129.91505611802876
    @test GasDispersion._monin_obukhov(1.2, ClassD) == Inf
    @test GasDispersion._monin_obukhov(1.2, ClassE) ≈ 129.91505611802876
    @test GasDispersion._monin_obukhov(1.2, ClassF) ≈ 26.818480014823884
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
        @test GasDispersion.crosswind_dispersion(1.2, class, DefaultSet()) ≈ cwind
        @test GasDispersion.vertical_dispersion(1.2, class, DefaultSet()) ≈ vert
    end

    # Puff dispersion
    knowns = [(ClassA, 0.21287234847910436, 0.6879188103871441),
              (ClassB, 0.16556738215041453, 0.6054498545741843),
              (ClassC, 0.11826241582172466, 0.3869880921561032),
              (ClassD, 0.07095744949303479, 0.17041904657983328),
              (ClassE, 0.047304966328689864, 0.11258170198247626),
              (ClassF, 0.023523465599668385, 0.05588182287353654)]
    @testset "Stability class $class" for (class,cwind,vert) in knowns
        @test GasDispersion.crosswind_dispersion(1.2, class, DefaultPuffSet()) ≈ cwind
        @test GasDispersion.downwind_dispersion(1.2, class, DefaultPuffSet()) ≈ cwind
        @test GasDispersion.vertical_dispersion(1.2, class, DefaultPuffSet()) ≈ vert
    end
end

@testset "Windspeed by powerlaw" begin

    u0, z0, p = 3.0, 1.0, 0.108
    a = SimpleAtmosphere(windspeed=u0, windspeed_height=z0, stability=ClassA)
    s = Scenario(Substance(:null,0,0,0,0,0,0,0,0,0,0,0),HorizontalJet(0,0,0,0,1.0,0,0,0),a)
    @test GasDispersion._windspeed(s) == GasDispersion._windspeed(a) ≈ u0
    @test GasDispersion._windspeed(s,10) == GasDispersion._windspeed(a,10) == GasDispersion._windspeed(u0,z0,10,ClassA,GasDispersion.DefaultWind)

    knowns = [(ClassA, 3.846991747968065),
              (ClassB, 3.882587524349958),
              (ClassC, 3.954770215669221),
              (ClassD, 4.160267486615665),
              (ClassE, 4.7876374417101975),
              (ClassF, 5.371817562105883)]

    @testset "Stability class $class" for (class, ans) in knowns
        @test  GasDispersion._windspeed(u0,z0,10,class,GasDispersion.DefaultWind) ≈ ans
    end

end

@testset "Windspeed by Monin-Obukhov length" begin
    ustar, zR, k = 3.0, 1.0, 0.35
    λ = GasDispersion._monin_obukhov(zR, ClassA)
    a(z) = (1-15*(z/λ))^0.25
    Ψ(a) = 2*log((1+a)/2) + log((1+a^2)/2) - 2*atan(a) + π/2

    @test GasDispersion._windspeed(10, ustar, zR, λ, ClassA, GasDispersion.BusingerWind) ≈ (ustar/k)*(log(10/zR) - Ψ(a(10)))
    @test GasDispersion._windspeed(10, ustar, zR, λ, ClassD, GasDispersion.BusingerWind) ≈ (ustar/k)*(log(10/zR))
    @test GasDispersion._windspeed(10, ustar, zR, λ, ClassE, GasDispersion.BusingerWind) ≈ (ustar/k)*(log(10/zR) + 4.7*(10/λ))
end

@testset "Plume rise" begin

    # test null case
    @test isa(GasDispersion.NoPlumeRise(),GasDispersion.PlumeRise)

    # test Briggs model of plume rise
    g = 9.80616 # m/s^2
    Tₐ = 288.15 # ambient temperature, K
    u = 4 # windspeed, m/s
    x = 50 # test point, m

    # Buoyant unstable plume
    # Fb = 50 case
    Tᵣ = 400 # K
    uⱼ = 72.93819699672669 # m/s
    xf = 565.0050541491846 # m
    Δhf = 100.71365158671999 # m
    sln = GasDispersion.plume_rise(1,uⱼ,Tᵣ,u,Tₐ,0,ClassA)
    @test isa(sln,GasDispersion.BuoyantPlume)
    @test sln ≈ GasDispersion.BuoyantPlume(50,xf,u,Δhf)
    @test GasDispersion.plume_rise(x, sln) ≈ 20.0
    @test GasDispersion.plume_rise(2xf, sln) ≈ Δhf

    # Fb = 60 case
    uⱼ = 87.52583639607202
    xf = 612.07897481385
    Δhf = 112.8895989623143
    sln = GasDispersion.plume_rise(1,uⱼ,Tᵣ,u,Tₐ,0,ClassA)
    @test isa(sln,GasDispersion.BuoyantPlume)
    @test sln ≈ GasDispersion.BuoyantPlume(60,xf,u,Δhf)

    # Momentum dominated unstable plume
    # Fb = 50 case
    Tᵣ = 325 # K
    uⱼ = 179.87751923862007 # m/s
    xf = 565.0050541491846 # m
    Fm = 7171.814541070672
    β = (1/3) + (u/uⱼ)
    Δhf = 3*(uⱼ/u) # m
    sln = GasDispersion.plume_rise(1,uⱼ,Tᵣ,u,Tₐ,0,ClassA)
    @test isa(sln,GasDispersion.MomentumPlume)
    @test sln ≈ GasDispersion.MomentumPlume(Fm,xf,β,u,0,Δhf,ClassA)
    @test GasDispersion.plume_rise(x, sln) ≈ 81.01824072514351
    @test GasDispersion.plume_rise(2xf, sln) ≈ Δhf

    # Buoyant stable plume
    # Fb = 50, class E
    Tᵣ = 400 # K
    uⱼ = 72.93819699672669 # m/s
    xf = 317.60677324769046 # m
    Δhf = 68.59722859012221 # m
    sln = GasDispersion.plume_rise(1,uⱼ,Tᵣ,u,Tₐ,0.02,ClassE)
    @test isa(sln,GasDispersion.BuoyantPlume)
    @test sln ≈ GasDispersion.BuoyantPlume(50,xf,u,Δhf)

    # Fb = 50, class F
    Tᵣ = 400 # K
    uⱼ = 72.93819699672669 # m/s
    xf = 240.0881533494489 # m
    Δhf = 56.92380039947288 # m
    sln = GasDispersion.plume_rise(1,uⱼ,Tᵣ,u,Tₐ,0.035,ClassF)
    @test isa(sln,GasDispersion.BuoyantPlume)
    @test sln ≈ GasDispersion.BuoyantPlume(50,xf,u,Δhf)

    # Momentum dominated stable plume
    # Fb = 67, class E
    Tᵣ = 325 # K
    uⱼ = 240.91531595745576 # m/s
    Fm = 12864.831225945449
    xf = 240.83782417699823 # m
    β = (1/3)+(u/uⱼ)
    s = 0.0006806288391462781
    Δhf = 19.04522445600697 # m
    sln = GasDispersion.plume_rise(1,uⱼ,Tᵣ,u,Tₐ,0.02,ClassE)
    @test isa(sln,GasDispersion.MomentumPlume)
    @test sln ≈ GasDispersion.MomentumPlume(Fm,xf,β,u,s,Δhf,ClassE)
    @test GasDispersion.plume_rise(x, sln) ≈ Δhf
    @test GasDispersion.plume_rise(2xf, sln) ≈ Δhf

    # Momentum dominated stable plume
    # Fb = 67, class F
    Tᵣ = 325 # K
    uⱼ = 240.91531595745576 # m/s
    Fm = 12864.831225945449
    xf = 182.0562825914961 # m
    β = (1/3)+(u/uⱼ)
    s = 0.0011911004685059867
    Δhf = 17.349211998937378 # m
    sln = GasDispersion.plume_rise(1,uⱼ,Tᵣ,u,Tₐ,0.035,ClassF)
    @test isa(sln,GasDispersion.MomentumPlume)
    @test sln ≈ GasDispersion.MomentumPlume(Fm,xf,β,u,s,Δhf,ClassF)

end

@testset "Equation Sets" begin
    
    include("ccps_tests.jl")
    include("tno_tests.jl")
    include("turner_tests.jl")
    include("isc3_tests.jl")

end