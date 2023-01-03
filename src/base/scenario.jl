# Substance type definition
struct Substance{D<:Union{Number,Function}, H<:Union{Number,Function}}
    name::Union{AbstractString,Symbol}
    ρ_g::D        # gas density, kg/m^3
    ρ_l::D        # liquid density, kg/m^3
    T_ref::Number # reference temperature for densities, K (default 15°C)
    P_ref::Number # reference pressure for densities, Pa (default 1atm)
    Tb::Number    # normal boiling point, K
    Δh_v::H       # specific enthalpy of vapourization, J/kg
    Cp_g::H       # specific heat capacity, gas, J/kg/K
    Cp_l::H       # specific heat capacity, liquid, J/kg/K
end
Substance(;name,gas_density,liquid_density,reference_temp=288.15,
reference_pressure=101325.0,boiling_temp,latent_heat,gas_heat_capacity,
liquid_heat_capacity) = Substance(name,gas_density,liquid_density,
reference_temp,reference_pressure,boiling_temp,latent_heat,gas_heat_capacity,
liquid_heat_capacity)

# Release type definition
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

# Stability classes
struct ClassA <: StabilityClass end
struct ClassB <: StabilityClass end
struct ClassC <: StabilityClass end
struct ClassD <: StabilityClass end
struct ClassE <: StabilityClass end
struct ClassF <: StabilityClass end

# Ambient atmosphere type definition
struct Ambient{T<:StabilityClass} <: Atmosphere
    pressure::Number          # Ambient pressure, Pa
    temperature::Number       # Ambient temperature, K
    density::Number           # Ambient density, kg/m^3
    windspeed::Number         # windspeed at windspeed height, m/s
    windspeed_height::Number  # reference height for windspeed, m
    stability::Type{T} # Pasquill-Gifford stability class
end
Ambient(; pressure=101325, temperature=298.15, density=1.225, windspeed=1.5,
        windspeed_height=10, stability=ClassF) = Ambient(pressure,temperature,
        density,windspeed,windspeed_height,stability)

# Scenario type definition
struct Scenario
    substance::Substance
    release::Release
    atmosphere::Atmosphere
end
Scenario(;substance, release, atmosphere=Ambient()) = Scenario(substance,
release, atmosphere)
