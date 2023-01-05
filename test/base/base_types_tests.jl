# I have shamelessly stolen this from the tests for show()
replstr(x, kv::Pair...) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80), kv...), MIME("text/plain"), x), x)

@testset "Base types tests" begin

@testset "Substance type" begin
    sub1 = Substance("test",8.90,490.0,298.15,120935.0368,100.,123456.,
     78910.,4.19)
    sub2 = Substance(name="test",gas_density=8.90,liquid_density=490.0,
     reference_temp=298.15,reference_pressure=120935.0368,boiling_temp=100.,
     latent_heat=123456.,gas_heat_capacity=78910.,liquid_heat_capacity=4.19)
    @test isa(sub1, Substance)
    @test sub1 ≈ sub2
    @test replstr(sub1) == "Substance: test \n"
end

@testset "Release type" begin
    rel1 = Release(1.0, 10.0, 0.25, 15.67, 2.0, 101325.0, 450., 0.67)
    rel2 = Release(mass_rate=1.0, duration=10.0, diameter=0.25,
     velocity=15.67, height=2.0, pressure=101325.0, temperature=450.,
     fraction_liquid=0.67)
    @test isa(rel1, Release)
    @test rel1 ≈ rel2
    @test replstr(rel1) == "Release conditions:\n    ṁ: 1.0 kg/s \n    Δt: 10.0 s \n    d: 0.25 m \n    u: 15.67 m/s \n    h: 2.0 m \n    P: 101325.0 Pa \n    T: 450.0 K \n    f_l: 0.67  \n"
end

@testset "DryAir type" begin
    atm1 = DryAir(100e3,273.15,287.05,2.0,5.0,ClassA)
    atm2 = DryAir(pressure=100e3,temperature=273.15,gas_constant=287.05,
     windspeed=2.0,windspeed_height=5.0,stability=ClassA)
    @test isa(atm1, DryAir)
    @test isa(atm1, Atmosphere)
    @test atm1 ≈ atm2
    @test replstr(atm1) == "Atmospheric conditions:\n    P: 100000.0 Pa \n    T: 273.15 K \n    Rs: 287.05 J/kg/K \n    u: 2.0 m/s \n    h: 5.0 m \n    stability: ClassA  \n"
end

@testset "Scenario type" begin
    sub = Substance("test",8.90,490.0,298.15,120935.0368,100.,123456., 78910.,4.19)
    rel = Release(1.0, 10.0, 0.25, 15.67, 2.0, 101325.0, 450., 0.67)
    atm = DryAir(100e3,273.15,287.05,2.0,5.0,ClassA)
    scn1 = Scenario(sub,rel,atm)
    scn2 = Scenario(substance=sub,release=rel,atmosphere=atm)
    @test isa(scn1, Scenario)
    @test scn1 ≈ scn2
    @test replstr(scn1) == replstr(sub) * replstr(rel) * replstr(atm)
end

end
