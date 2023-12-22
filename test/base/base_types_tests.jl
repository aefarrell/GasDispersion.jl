# I have shamelessly stolen this from the tests for show()
replstr(x, kv::Pair...) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80), kv...), MIME("text/plain"), x), x)

@testset "Base types tests" begin

@testset "Substance type" begin
    sub1 = Substance("test",2.0,1.0,8.90,490,298.15,120935.0368,1.3,100,123456,78910,4.19)
    sub2 = Substance(name="test",molar_weight=2.0,vapor_pressure=1.0,gas_density=8.90,liquid_density=490,
     reference_temp=298.15,reference_pressure=120935.0368,k=1.3,boiling_temp=100,
     latent_heat=123456,gas_heat_capacity=78910,liquid_heat_capacity=4.19)
    @test isa(sub1, Substance)
    @test sub1 ≈ sub2
    @test replstr(sub1) == "Substance: test \n"
end

@testset "Release type" begin
    rel1 = HorizontalJet(1, 10, 0.25, 15.67, 2, 101325, 450, 0.67)
    rel2 = HorizontalJet(mass_rate=1, duration=10, diameter=0.25,
     velocity=15.67, height=2, pressure=101325, temperature=450,
     fraction_liquid=0.67)
    @test isa(rel1, Release)
    @test rel1 ≈ rel2
    @test replstr(rel1) == "HorizontalJet release:\n    ṁ: 1.0 kg/s \n    Δt: 10.0 s \n    d: 0.25 m \n    u: 15.67 m/s \n    h: 2.0 m \n    P: 101325.0 Pa \n    T: 450.0 K \n    f_l: 0.67  \n"
end

@testset "SimpleAtmosphere type" begin
    atm1 = SimpleAtmosphere(100e3,273.15,2,5,0,ClassA)
    atm2 = SimpleAtmosphere(pressure=100e3,temperature=273.15,windspeed=2,
                            windspeed_height=5,stability=ClassA)
    @test isa(atm1, SimpleAtmosphere)
    @test isa(atm1, Atmosphere)
    @test atm1 ≈ atm2
    @test isnothing(GasDispersion._lapse_rate(atm1))
    @test GasDispersion._lapse_rate(SimpleAtmosphere(stability=ClassE)) == 0.020
    @test GasDispersion._lapse_rate(SimpleAtmosphere(stability=ClassF)) == 0.035
    @test GasDispersion._rel_humidity(atm1) == 0.0
    @test replstr(atm1) == "SimpleAtmosphere atmosphere:\n    P: 100000.0 Pa \n    T: 273.15 K \n    u: 2.0 m/s \n    h: 5.0 m \n    rh: 0.0 % \n    stability: ClassA  \n"
end

@testset "Scenario type" begin
    sub = Substance("test",2.0,1.0,8.90,490,298.15,120935.0368,1.3,100,123456,78910,4.19)
    rel = HorizontalJet(1, 10, 0.25, 15.67, 2, 101325, 450, 0.67)
    atm = SimpleAtmosphere(100e3,273.15,2,5,0,ClassA)
    scn1 = Scenario(sub,rel,atm)
    scn2 = Scenario(substance=sub,release=rel,atmosphere=atm)
    @test isa(scn1, Scenario)
    @test scn1 ≈ scn2
    @test replstr(scn1) == replstr(sub) * replstr(rel) * replstr(atm)
end

@testset "Property getters" begin
    sub = Substance("test",2.0,1.0,8.90,490,298.15,120935.0368,1.3,100,123456,78910,4.19)
    rel = HorizontalJet(1, 10, 0.25, 15.67, 2, 101325, 450, 0.67)
    atm = SimpleAtmosphere(100e3,273.15,2,5,0,ClassA)
    scn = Scenario(sub,rel,atm)

    @test GasDispersion._atmosphere_temperature(scn) == GasDispersion._temperature(atm) == 273.15
    @test GasDispersion._release_temperature(scn) == GasDispersion._temperature(rel) == 450.0
    @test GasDispersion._atmosphere_pressure(scn) == GasDispersion._pressure(atm) == 100e3
    @test GasDispersion._release_pressure(scn) == GasDispersion._pressure(rel) == 101325.0

    @test GasDispersion._mass_rate(scn) == GasDispersion._mass_rate(rel) == 1.0
    @test GasDispersion._duration(scn) == GasDispersion._duration(rel) == 10.0
    @test GasDispersion._release_mass(scn) == GasDispersion._mass(rel) == 10.0
    @test GasDispersion._release_diameter(scn) == GasDispersion._diameter(rel) == 0.25
    @test GasDispersion._release_area(scn) == GasDispersion._area(rel) ≈ (π/4)*0.25^2
    @test GasDispersion._release_velocity(scn) == GasDispersion._velocity(rel) == 15.67
    @test GasDispersion._release_flowrate(scn) == GasDispersion._flowrate(rel) ≈ (π/4)*(0.25^2)*15.67
    @test GasDispersion._release_height(scn) == GasDispersion._height(rel) == 2.0
    @test GasDispersion._release_liquid_fraction(scn) == GasDispersion._liquid_fraction(rel) == 0.67

    @test GasDispersion._windspeed(scn) == GasDispersion._windspeed(atm) == 2.0
    @test GasDispersion._windspeed_height(scn) == GasDispersion._windspeed_height(atm) == 5.0
    @test GasDispersion._stability(scn) == GasDispersion._stability(atm) == ClassA

end

@testset "Density functions" begin
    sub1 = Substance("test",2.0,1.0,8.90,490,298.15,120935.0368,1.3,100,123456,78910,4.19)
    sub2 = Substance("test",2.0,1.0,(x,y)->y*x^2,(x,y)->y*x^3,2,3,1.3,100,123456., 78910.,4.19)
    rel = HorizontalJet(1.0, 10.0, 0.25, 15.67, 2.0, 101325.0, 450., 0.5)
    atm = SimpleAtmosphere(100e3,273.15,2.0,5.0,0.0,ClassA)
    scn1 = Scenario(sub1,rel,atm)
    scn2 = Scenario(sub2,rel,atm)

    @test GasDispersion._liquid_density(sub1) ≈ 490.0
    @test GasDispersion._liquid_density(sub2) ≈ 24
    @test GasDispersion._gas_density(sub1) ≈ 8.90
    @test GasDispersion._gas_density(sub2) ≈ 12

    @test GasDispersion._density(sub1, 0.5, 298.15, 120935.0368) ≈ 2/(1/8.9 + 1/490.)
    @test GasDispersion._density(sub2, 0.5, 2, 3) ≈ 48/3

    @test GasDispersion._release_density(scn1) == GasDispersion._density(sub1,0.5,450.0,101325.0)
    @test GasDispersion._release_density(scn2) == GasDispersion._density(sub2,0.5,450.0,101325.0)

    @test GasDispersion._density(atm, 100, 200) ≈ 2/(8.31446261815324/0.028960)
    @test GasDispersion._atmosphere_density(scn1) == GasDispersion._density(atm)

end

end
