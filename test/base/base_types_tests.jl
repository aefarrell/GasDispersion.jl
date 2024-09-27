# I have shamelessly stolen this from the tests for show()
replstr(x, kv::Pair...) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80), kv...), MIME("text/plain"), x), x)

@testset "Base types tests" begin

@testset "Substance type" begin

    sub = Substance("test",2.0,1.0,8.90,490,298.15,120935.0368,1.3,100,123456,78910,4.19)
    @test sub isa Substance
    @test replstr(sub) == "Substance: test \n    MW: 2.0 kg/mol \n    P_v: 1.0 Pa \n    ρ_g: 8.9 kg/m^3 \n    ρ_l: 490 kg/m^3 \n    T_ref: 298.15 K \n    P_ref: 120935.0368 Pa \n    k: 1.3  \n    T_b: 100.0 K \n    Δh_v: 123456 J/kg \n    Cp_g: 78910 J/kg/K \n    Cp_l: 4.19 J/kg/K \n"

    R = GasDispersion.R_GAS_CONST

    # Antoine equation, propane, *Properties of Gases and Liquids, 5th ed*
    # constants converted from log_10(Pv [bar]) to ln(Pv [Pa])
    pv_ant = GasDispersion.Antoine(log(10)*(3.92828+5),log(10)*803.9970,247.040-273.15)
    @test pv_ant(247.76) ≈ 199964.60671380349

    # propane from Perry's 8th edition, DIPPR correlations
    pv = GasDispersion.DIPPRVaporPressure(59.078,-3_492.6,-6.0669,1.0919e-5,2)
    ρl = GasDispersion.DIPPRLiquidDensity(0.044096,369.83,1.3757,0.27453,369.83,0.29359)
    ρg = 0.044096*101325/(R*288.15)
    Δhv = GasDispersion.DIPPRLatentHeat(0.044096,369.83,2.9209e7,0.78237,-0.77319,0.39246,0)
    cp_ig = GasDispersion.DIPPRIdealGasHeatCapacity(0.044096,369.83,0.5192e5,1.9245e5,1.6265e3,1.168e5,723.6)
    cp_l = GasDispersion.DIPPRLiquidHeatCapacity(GasDispersion.Eq2,0.044096,369.83,62.983,113_630,633.21,-873.46,0)

    propane = Substance(name="propane",
                        molar_weight=0.044096,
                        vapor_pressure=pv,
                        gas_density=nothing,
                        liquid_density=ρl,
                        reference_temp=288.15,
                        reference_pressure=101325,
                        k=1.3,
                        boiling_temp=231.02,
                        latent_heat=Δhv,
                        gas_heat_capacity=cp_ig,
                        liquid_heat_capacity=cp_l)
    # @test isa(propane, Substance)
    # @test replstr(propane) == "Substance: propane \n    MW: 0.044096 kg/mol \n    P_v: GasDispersion.DIPPRVaporPressure{Float64}(59.078, -3492.6, -6.0669, 1.0919e-5, 2.0) Pa \n    ρ_g: 1.864931992847327 kg/m^3 \n    ρ_l: GasDispersion.DIPPRLiquidDensity{Float64}(0.044096, 369.83, 1.3757, 0.27453, 369.83, 0.29359) kg/m^3 \n    T_ref: 288.15 K \n    P_ref: 101325.0 Pa \n    k: 1.3  \n    T_b: 231.02 K \n    Δh_v: GasDispersion.DIPPRLatentHeat{Float64}(0.044096, 369.83, 2.9209e7, 0.78237, -0.77319, 0.39246, 0.0) J/kg \n    Cp_g: GasDispersion.DIPPRIdealGasHeatCapacity{Float64}(0.044096, 369.83, 51920.0, 192450.0, 1626.5, 116800.0, 723.6) J/kg/K \n    Cp_l: GasDispersion.DIPPRLiquidHeatCapacity{GasDispersion.Eq2, Float64}(GasDispersion.Eq2, 0.044096, 369.83, 62.983, 113630.0, 633.21, -873.46, 0.0) J/kg/K \n"

    @test GasDispersion._MW(propane) ≈ 0.044096
    @test GasDispersion._boiling_temperature(propane) ≈ 231.02
    @test GasDispersion._vapor_pressure(propane,288.15) == pv(288.15) ≈ 732387.3305651173
    @test GasDispersion._liquid_density(propane) == GasDispersion._liquid_density(propane,288.15,101325) == ρl(288.15,101325) ≈ 506.61518872042035
    @test GasDispersion._gas_density(propane) == GasDispersion._gas_density(propane,288.15,101325) ≈ ρg
    @test GasDispersion._density(propane,0.5,288.15,101325) ≈ 2/(1/ρg + 1/506.61518872042035)
    @test GasDispersion._latent_heat(propane) == GasDispersion._latent_heat(propane,231.02) == Δhv(231.02) ≈ 425139.0953430851
    @test GasDispersion._cp_gas(propane) == GasDispersion._cp_gas(propane,288.15) == cp_ig(288.15) ≈ 1618.86346283786
    @test GasDispersion._cp_liquid(propane) == GasDispersion._cp_liquid(propane,288.15) == cp_l(288.15) ≈ 2626.002694859999

    Base.isapprox(a::GasDispersion.Antoine, b::GasDispersion.Antoine) = all([
        getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
        if typeof(getproperty(a,k))<:Number ])

    test_pv = Substance(name="propane",
                        molar_weight=0.044096,
                        vapor_pressure=nothing,
                        gas_density=nothing,
                        liquid_density=ρl,
                        reference_temp=288.15,
                        reference_pressure=101325,
                        k=1.142,
                        boiling_temp=231.02,
                        latent_heat=Δhv,
                        gas_heat_capacity=cp_ig,
                        liquid_heat_capacity=cp_l)
    @test test_pv.P_v ≈ GasDispersion.Antoine(9.759924888223345, 2254.7378476773574, 0.0)
    
    test_pv2 = Substance(name="propane",
                         molar_weight=0.044096,
                         vapor_pressure=nothing,
                         gas_density=nothing,
                         liquid_density=ρl,
                         reference_temp=288.15,
                         reference_pressure=101325,
                         k=1.142,
                         boiling_temp=231.02,
                         latent_heat=Δhv(288.15),
                         gas_heat_capacity=cp_ig,
                         liquid_heat_capacity=cp_l)
    @test test_pv2.P_v ≈ GasDispersion.Antoine(8.086226333190652, 1868.0800074937044, 0.0)

    # adding this to be tediously completionist
    @test GasDispersion.DIPPRLiquidDensity(0.044096,369.83,1,2,3,4) == GasDispersion.DIPPRLiquidDensity(0.044096,369.83,1.0,2.0,3.0,4.0)

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
    # check that the air and water correlations are correct
    # i.e. match known values from Perry's
    @test GasDispersion._MW(GasDispersion.DRYAIR) ≈ 0.028960
    @test GasDispersion._vapor_pressure(GasDispersion.DRYAIR,132.45) ≈ 3.794896199565284e6 # Table 2-8
    @test GasDispersion._gas_density(GasDispersion.DRYAIR) ≈ 1.2247920562603996 # Ideal Gas Law
    @test GasDispersion._liquid_density(GasDispersion.DRYAIR, 132.45, 101325)/28.96 ≈ 10.83417498971309 # Table 2-32
    @test GasDispersion._latent_heat(GasDispersion.DRYAIR, 59.15)*28.96 ≈ 6.759007660367378e6 # Table 2-150
    @test GasDispersion._cp_gas(GasDispersion.DRYAIR, 50)*28.96 ≈ 28958.0 # Table 2-156
    @test GasDispersion._cp_liquid(GasDispersion.DRYAIR, 75)*28.96 ≈ 53065.0 # Table 2-153

    @test GasDispersion._MW(GasDispersion.WATER) ≈ 0.018015
    @test GasDispersion._vapor_pressure(GasDispersion.WATER, 273.16) ≈  610.5626315257215 # Table 2-8
    @test GasDispersion._gas_density(GasDispersion.WATER) ≈  0.7619001689755213 # Ideal Gas Law
    @test GasDispersion._liquid_density(GasDispersion.WATER, 273.16, 101325)/18.015 ≈  53.07493666634459 # Hand calc, equation given in Table 2-32
    @test GasDispersion._latent_heat(GasDispersion.WATER, 273.16)*18.015 ≈ 4.473232632794751e7 # Table 2-150
    @test GasDispersion._cp_gas(GasDispersion.WATER, 100)*18.015 ≈ 33363.000341254825 # Table 2-156
    @test GasDispersion._cp_liquid(GasDispersion.WATER, 273.16)*18.015 ≈ 76150.12956433237 # Table 2-153

    atm1 = SimpleAtmosphere(100e3,273.15,2,5,10,ClassA)
    atm2 = SimpleAtmosphere(pressure=100e3,temperature=273.15,windspeed=2,
                            windspeed_height=5,rel_humidity=10,stability=ClassA)
    @test isa(atm1, SimpleAtmosphere)
    @test isa(atm1, Atmosphere)
    @test atm1 ≈ atm2
    @test replstr(atm1) == "SimpleAtmosphere atmosphere:\n    P: 100000.0 Pa \n    T: 273.15 K \n    u: 2.0 m/s \n    h: 5.0 m \n    rh: 10.0 % \n    stability: ClassA  \n"

    @test isnothing(GasDispersion._lapse_rate(atm1))
    @test GasDispersion._lapse_rate(SimpleAtmosphere(stability=ClassE)) == 0.020
    @test GasDispersion._lapse_rate(SimpleAtmosphere(stability=ClassF)) == 0.035
    @test GasDispersion._rel_humidity(atm1) ≈ 10.0
    @test GasDispersion._surface_roughness(atm1) == 1.0
    @test GasDispersion._density(atm1) ≈ 1.2748615244572732
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

    @test GasDispersion._release_density(scn1) == GasDispersion._density(sub1,0.5,450.0,101325.0)
    @test GasDispersion._release_density(scn2) == GasDispersion._density(sub2,0.5,450.0,101325.0)

    @test GasDispersion._density(atm, 100, 200) ≈ 2/(8.31446261815324/0.028960)
    @test GasDispersion._atmosphere_density(scn1) == GasDispersion._density(atm)

end

end
