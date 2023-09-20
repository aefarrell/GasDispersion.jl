# units for the show function
const units = Dict{Symbol,String}([
    :ṁ => "kg/s",
    :Δt => "s",
    :d => "m",
    :u => "m/s",
    :h => "m",
    :ρ => "kg/m^3",
    :P => "Pa",
    :T => "K",
    :Rs => "J/kg/K",
    :stability => "",
    :f => ""
])

# Substance type definition
struct Substance{D_G,D_L,H,CP_G,CP_L}
    name::Union{AbstractString,Symbol}
    ρ_g::D_G      # gas density, kg/m^3
    ρ_l::D_L      # liquid density, kg/m^3
    T_ref::Number # reference temperature for densities, K (default 15°C)
    P_ref::Number # reference pressure for densities, Pa (default 1atm)
    T_b::Number   # normal boiling point, K
    Δh_v::H       # specific enthalpy of vapourization, J/kg
    Cp_g::CP_G    # specific heat capacity, gas, J/kg/K
    Cp_l::CP_L    # specific heat capacity, liquid, J/kg/K
end
Substance(;name,gas_density,liquid_density,reference_temp=288.15,
reference_pressure=101325.0,boiling_temp,latent_heat,gas_heat_capacity,
liquid_heat_capacity) = Substance(name,gas_density,liquid_density,
reference_temp,reference_pressure,boiling_temp,latent_heat,gas_heat_capacity,
liquid_heat_capacity)

Base.isapprox(a::Substance, b::Substance) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

function Base.show(io::IO, mime::MIME"text/plain", s::Substance)
    print(io, "Substance: $(s.name) \n")
end

# Release type definition
struct Release
    ṁ::Number       # mass emission rate, kg/s
    Δt::Number        # release duration, s
    d::Number        # release diameter, m
    u::Number        # release velocity, m/s
    h::Number          # release height, m
    P::Number        # release pressure, Pa
    T::Number     # release temperature, K
    f_l::Number #  mass fraction liquid, unitless
end
Release(; mass_rate, duration, diameter, velocity, height, pressure,
    temperature, fraction_liquid) = Release(mass_rate, duration, diameter,
    velocity, height, pressure, temperature, fraction_liquid)

Base.isapprox(a::Release, b::Release) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

function Base.show(io::IO, mime::MIME"text/plain", r::Release)
    print(io, "Release conditions:\n")
    for key in fieldnames(typeof(r))
        val =  getproperty(r, key)
        var =  Symbol(split(string(key),"_")[1])
        unit = units[var]
        print(io, "    $key: $val $unit \n")
    end
end


# Stability classes
struct ClassA <: StabilityClass end
struct ClassB <: StabilityClass end
struct ClassC <: StabilityClass end
struct ClassD <: StabilityClass end
struct ClassE <: StabilityClass end
struct ClassF <: StabilityClass end

# DryAir atmosphere type definition
struct DryAir{T<:StabilityClass} <: Atmosphere
    P::Number  # atmospheric pressure, Pa
    T::Number  # atmospheric temperature, K
    Rs::Number # specific gas constant for dry air, 287.0500676 J/kg/K
    u::Number  # windspeed at windspeed height, m/s
    h::Number  # reference height for windspeed, m
    stability::Type{T} # Pasquill-Gifford stability class
end
DryAir(; pressure=101325,temperature=298.15,gas_constant=287.0500676,
        windspeed=1.5,windspeed_height=10, stability=ClassF) = DryAir(pressure,
        temperature,gas_constant,windspeed,windspeed_height,stability)

Base.isapprox(a::Atmosphere, b::Atmosphere) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

function Base.show(io::IO, mime::MIME"text/plain", a::Atmosphere)
    print(io, "Atmospheric conditions:\n")
    for key in fieldnames(typeof(a))
        val =  getproperty(a, key)
        unit = units[key]
        print(io, "    $key: $val $unit \n")
    end
end

# Scenario type definition
struct Scenario
    substance::Substance
    release::Release
    atmosphere::Atmosphere
end
Scenario(;substance, release, atmosphere=DryAir()) = Scenario(substance,
release, atmosphere)

Base.isapprox(a::Scenario, b::Scenario) = all( [ a.substance≈b.substance,
                                                 a.release≈b.release,
                                                 a.atmosphere≈b.atmosphere ])

function Base.show(io::IO, mime::MIME"text/plain", s::Scenario)
 show(io,mime,s.substance)
 show(io,mime,s.release)
 show(io,mime,s.atmosphere)
end
