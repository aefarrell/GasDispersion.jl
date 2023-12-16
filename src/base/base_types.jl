# units for the show function
units = Dict{Symbol,String}([
    :ṁ => "kg/s",
    :Δt => "s",
    :d => "m",
    :u => "m/s",
    :h => "m",
    :ρ => "kg/m^3",
    :P => "Pa",
    :T => "K",
    :Rs => "J/kg/K",
    :rh => "%",
    :stability => "",
    :f => ""
])

# Substance type definition
struct Substance{N<:Union{AbstractString,Symbol},D_G,D_L,F<:Number,H,CP_G,CP_L}
    name::N
    ρ_g::D_G    # gas density, kg/m^3
    ρ_l::D_L    # liquid density, kg/m^3
    T_ref::F    # reference temperature for densities, K (default 15°C)
    P_ref::F    # reference pressure for densities, Pa (default 1atm)
    T_b::F      # normal boiling point, K
    Δh_v::H     # specific enthalpy of vapourization, J/kg
    Cp_g::CP_G  # specific heat capacity, gas, J/kg/K
    Cp_l::CP_L  # specific heat capacity, liquid, J/kg/K
end
Substance(name,ρ_g,ρ_l,T_ref,P_ref,T_b,Δh_v,Cp_g,Cp_l) = Substance(name,ρ_g,ρ_l,
    promote(T_ref,P_ref,T_b)...,Δh_v,Cp_g,Cp_l) 
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
struct HorizontalJet{F <: Number} <: Release
    ṁ::F   # mass emission rate, kg/s
    Δt::F  # release duration, s
    d::F   # release diameter, m
    u::F   # release velocity, m/s
    h::F   # release height, m
    P::F   # release pressure, Pa
    T::F   # release temperature, K
    f_l::F #  mass fraction liquid, unitless
end
HorizontalJet(ṁ,Δt,d,u,h,P,T,f_l) = HorizontalJet(promote(ṁ,Δt,d,u,h,P,T,f_l)...)
HorizontalJet(; mass_rate, duration, diameter, velocity, height, pressure,
    temperature, fraction_liquid) = HorizontalJet(mass_rate, duration, diameter,
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

# SimpleAtmosphere atmosphere type definition
struct SimpleAtmosphere{F<:Number,S<:StabilityClass} <: Atmosphere
    P::F  # atmospheric pressure, Pa
    T::F # atmospheric temperature, K
    Rs::F # specific gas constant for dry air, 287.0500676 J/kg/K
    u::F  # windspeed at windspeed height, m/s
    h::F  # reference height for windspeed, m
    rh::F      # relative humidity, %
    stability::Type{S} # Pasquill-Gifford stability class
end
SimpleAtmosphere(P,T,Rs,u,h,rh,stability) = SimpleAtmosphere(promote(P,T,Rs,u,h,rh,)...,stability)
SimpleAtmosphere(; pressure=101325,temperature=298.15,gas_constant=287.0500676,
        windspeed=1.5,windspeed_height=10,rel_humidity=0.0,stability=ClassF) = SimpleAtmosphere(pressure,
        temperature,gas_constant,windspeed,windspeed_height,rel_humidity,stability)

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
struct Scenario{S<:Substance,R<:Release,A<:Atmosphere}
    substance::S
    release::R
    atmosphere::A
end
Scenario(;substance, release, atmosphere=SimpleAtmosphere()) = Scenario(substance,
release, atmosphere)

Base.isapprox(a::Scenario, b::Scenario) = all( [ a.substance≈b.substance,
                                                 a.release≈b.release,
                                                 a.atmosphere≈b.atmosphere ])

function Base.show(io::IO, mime::MIME"text/plain", s::Scenario)
 show(io,mime,s.substance)
 show(io,mime,s.release)
 show(io,mime,s.atmosphere)
end

# Default equation set
struct DefaultSet <: EquationSet end