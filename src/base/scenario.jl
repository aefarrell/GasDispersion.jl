struct Release
    mass_rate::Number       # mass emission rate, kg/s
    duration::Number        # release duration, s
    diameter::Number        # release diameter, m
    velocity::Number        # release velocity, m/s
    height::Number          # release height, m
    pressure::Number        # release pressure, Pa
    temperature::Number     # release temperature, K
    fraction_liquid::Number #  mass fraction liquid, unitless
end
Release(; mass_rate, duration, diameter, velocity, height, pressure,
    temperature, fraction_liquid) = Release(mass_rate, duration, diameter,
    velocity, height, pressure, temperature, fraction_liquid)

struct Substance{D<:Union{Number,Function}, H<:Union{Number,Function}}
    name::Union{AbstractString,Symbol}
    ρ_g::D      # gas density, kg/m^3
    ρ_l::D      # liquid density, kg/m^3
    Tb::Number  # normal boiling point, K
    Δh_v::H     # specific enthalpy of vapourization, J/kg
    Cp_g::H     # specific heat capacity, gas, J/kg/K
    Cp_l::H     # specific heat capacity, liquid, J/kg/K
end
Substance(;name, gas_density, liquid_density, boiling_temp, latent_heat,
gas_heat_capacity, liquid_heat_capacity) = Substance(name, gas_density,
liquid_density, boiling_temp, latent_heat, gas_heat_capacity, liquid_heat_capacity)

struct Scenario
    substance::Substance
    release::Release
    atmosphere::Atmosphere
end
Scenario(;substance, release, atmosphere=Ambient()) = Scenario(substance,
release, atmosphere)

struct Ambient <: Atmosphere
    pressure::Number
    temperature::Number
    density::Number
    windspeed::Number
    stability::Symbol
end
Ambient(; pressure=101325, temperature=298.15, density=1.225, windspeed=1.5,
          stability=:F) = Ambient(pressure,temperature,density,windspeed,stability)
