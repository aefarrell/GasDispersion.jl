@testset "ClapeyronSubstance type" begin

    ethylene = PR("ethylene")
    sub = Substance(ethylene, reference_temp=169)
    @test sub isa AbstractSubstance

    @test sub.name == ethylene.components[1]

    T,P = 169, 101325
    MW = Clapeyron.mw(ethylene)[1]
    @test GasDispersion._MW(sub) == MW/1000
    @test GasDispersion._vapor_pressure(sub, T) ≈ Clapeyron.saturation_pressure(ethylene, T)[1]
    @test GasDispersion._cp_gas(sub) ≈ GasDispersion._cp_gas(sub,T) ≈ Clapeyron.isobaric_heat_capacity(ethylene, P, T; phase=:gas)*1000/MW
    @test GasDispersion._cp_liquid(sub) ≈ GasDispersion._cp_liquid(sub,T) ≈ Clapeyron.isobaric_heat_capacity(ethylene, P, T; phase=:liquid)*1000/MW
    @test GasDispersion._boiling_temperature(sub) ≈ Clapeyron.saturation_temperature(ethylene, 101325)[1]
    @test GasDispersion._latent_heat(sub,T) ≈ Clapeyron.enthalpy_vap(ethylene, T)*1000/MW
    @test GasDispersion._latent_heat(sub) ≈ Clapeyron.enthalpy_vap(ethylene, sub.T_b)*1000/MW

    # density functions
    @test GasDispersion._liquid_density(sub) ≈ GasDispersion._liquid_density(sub, T, P) ≈ Clapeyron.mass_density(ethylene,P,T; phase=:liquid)
    @test GasDispersion._gas_density(sub) ≈ GasDispersion._gas_density(sub, T, P) ≈ Clapeyron.mass_density(ethylene,P,T; phase=:gas)

end