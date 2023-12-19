# I have shamelessly stolen this from the tests for show()
replstr(x, kv::Pair...) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80), kv...), MIME("text/plain"), x), x)

@testset "Base types tests" begin

@testset "Substance type" begin
    sub1 = Substance("test",1.0,8.90,490,298.15,120935.0368,1.3,100,123456,78910,4.19)
    sub2 = Substance(name="test",vapor_pressure=1.0,gas_density=8.90,liquid_density=490,
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
    @test replstr(rel1) == "HorizontalJet{Float64} release:\n    ṁ: 1.0 kg/s \n    Δt: 10.0 s \n    d: 0.25 m \n    u: 15.67 m/s \n    h: 2.0 m \n    P: 101325.0 Pa \n    T: 450.0 K \n    f_l: 0.67  \n"
end

@testset "SimpleAtmosphere type" begin
    atm1 = SimpleAtmosphere(100e3,273.15,287.05,2,5,0,ClassA)
    atm2 = SimpleAtmosphere(pressure=100e3,temperature=273.15,gas_constant=287.05,
     windspeed=2,windspeed_height=5,stability=ClassA)
    @test isa(atm1, SimpleAtmosphere)
    @test isa(atm1, Atmosphere)
    @test atm1 ≈ atm2
    @test replstr(atm1) == "SimpleAtmosphere{Float64, ClassA} atmosphere:\n    P: 100000.0 Pa \n    T: 273.15 K \n    Rs: 287.05 J/kg/K \n    u: 2.0 m/s \n    h: 5.0 m \n    rh: 0.0 % \n    stability: ClassA  \n"
end

@testset "Scenario type" begin
    sub = Substance("test",1.0,8.90,490,298.15,120935.0368,1.3,100,123456,78910,4.19)
    rel = HorizontalJet(1, 10, 0.25, 15.67, 2, 101325, 450, 0.67)
    atm = SimpleAtmosphere(100e3,273.15,287.05,2,5,0,ClassA)
    scn1 = Scenario(sub,rel,atm)
    scn2 = Scenario(substance=sub,release=rel,atmosphere=atm)
    @test isa(scn1, Scenario)
    @test scn1 ≈ scn2
    @test replstr(scn1) == replstr(sub) * replstr(rel) * replstr(atm)
end

end
