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
    :f => "",
    :k => ""
])

struct Antoine{F<:Number}
    A::F
    B::F
    C::F
end

# Substance type definition
struct Substance{N<:Union{AbstractString,Symbol},MOL_WT,VAP,D_G,D_L,F<:Number,H,CP_G,CP_L}
    name::N
    MW::MOL_WT  # molar weight, kg/mol
    P_v::VAP    # vapour pressure, Pa
    ρ_g::D_G    # gas density, kg/m^3
    ρ_l::D_L    # liquid density, kg/m^3
    T_ref::F    # reference temperature for densities, K (default 15°C)
    P_ref::F    # reference pressure for densities, Pa (default 1atm)
    k::F        # heat capacity ratio Cp/Cv, unitless (default 1.4)
    T_b::F      # normal boiling point, K
    Δh_v::H     # specific enthalpy of vapourization, J/kg
    Cp_g::CP_G  # specific heat capacity, gas, J/kg/K
    Cp_l::CP_L  # specific heat capacity, liquid, J/kg/K
end

function Substance(name,molar_weight,vapor_pressure,gas_density,liquid_density,
                   reference_temp,reference_pressure,k,boiling_temp,latent_heat,
                   gas_heat_capacity,liquid_heat_capacity)

    R = 8.31446261815324
    
    # initialize with ideal gas
    if isnothing(gas_density)
        gas_density = reference_pressure*molar_weight/(R*reference_temp)
    end

    # if no vapour pressure correlation given, use Clapeyron equation
    if isnothing(vapor_pressure)
        B = latent_heat*molar_weight/R
        A = B/(boiling_temp)
        C = 0.0
        vapor_pressure = Antoine(A,B,C)
    end

    return Substance(name,molar_weight,vapor_pressure,gas_density,liquid_density,
    promote(reference_temp,reference_pressure,k,boiling_temp)...,latent_heat,gas_heat_capacity,
    liquid_heat_capacity)
end

Substance(;name,molar_weight,vapor_pressure=nothing,gas_density=nothing,liquid_density,reference_temp=288.15,
           reference_pressure=101325.0,k=1.4,boiling_temp,latent_heat,gas_heat_capacity,
           liquid_heat_capacity) = Substance(name,molar_weight,vapor_pressure,gas_density,liquid_density,
           reference_temp,reference_pressure,k,boiling_temp,latent_heat,gas_heat_capacity,
           liquid_heat_capacity)

Base.isapprox(a::Substance, b::Substance) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

function Base.show(io::IO, mime::MIME"text/plain", s::Substance)
    print(io, "Substance: $(s.name) \n")
end

# Substance property getters
_MW(s::Substance) = s.MW
_vapor_pressure(s::Substance{<:Any,<:Any,<:Antoine,<:Any,<:Any,<:Any,<:Any,<:Any}, T) = s.P_ref*exp(s.P_v.A - s.P_v.B/(T + s.P_v.C))
_vapor_pressure(s::Substance{<:Any,<:Any,<:Function,<:Any,<:Any,<:Any,<:Any,<:Any}, T) = s.P_ref*s.P_v(T)
_cp_gas(s::Substance{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Number,<:Any}) = s.Cp_g
_cp_liquid(s::Substance{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Number}) = s.Cp_l
_boiling_temperature(s::Substance) = s.T_b
_latent_heat(s::Substance{<:Any,<:Any,<:Any,<:Any,<:Any,<:Number,<:Any,<:Any}) = s.Δh_v

# density functions
_liquid_density(s::Substance) = _liquid_density(s, s.T_ref, s.P_ref)
_liquid_density(s::Substance{<:Any,<:Any,<:Any,<:Any,<:Number,<:Any,<:Any,<:Any,<:Any}, T::Number, P::Number) = s.ρ_l
_liquid_density(s::Substance{<:Any,<:Any,<:Any,<:Any,<:Function,<:Any,<:Any,<:Any,<:Any}, T::Number, P::Number) = s.ρ_l(T,P)

_gas_density(s::Substance) = _gas_density(s, s.T_ref, s.P_ref)
_gas_density(s::Substance{<:Any,<:Any,<:Any,<:Number,<:Any,<:Any,<:Any,<:Any,<:Any}, T::Number, P::Number) = s.ρ_g*(s.T_ref/T)*(P/s.P_ref)
_gas_density(s::Substance{<:Any,<:Any,<:Any,<:Function,<:Any,<:Any,<:Any,<:Any,<:Any}, T::Number, P::Number) = s.ρ_g(T,P)

function _density(s::Substance, f_l, T, P)
    f_g = 1 - f_l
    ρ_l = _liquid_density(s,T,P)
    ρ_g = _gas_density(s,T,P)

    return 1/(f_l/ρ_l + f_g/ρ_g)
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

struct VerticalJet{F <: Number} <: Release
    ṁ::F   # mass emission rate, kg/s
    Δt::F  # release duration, s
    d::F   # release diameter, m
    u::F   # release velocity, m/s
    h::F   # release height, m
    P::F   # release pressure, Pa
    T::F   # release temperature, K
    f_l::F #  mass fraction liquid, unitless
end
VerticalJet(ṁ,Δt,d,u,h,P,T,f_l) = VerticalJet(promote(ṁ,Δt,d,u,h,P,T,f_l)...)
VerticalJet(; mass_rate, duration, diameter, velocity, height, pressure,
    temperature, fraction_liquid) = VerticalJet(mass_rate, duration, diameter,
    velocity, height, pressure, temperature, fraction_liquid)

Base.isapprox(a::Release, b::Release) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

function Base.show(io::IO, mime::MIME"text/plain", r::Release)
    r_type = split(string(typeof(r)),"{")[1]
    print(io, "$r_type release:\n")
    for key in fieldnames(typeof(r))
        val =  getproperty(r, key)
        var =  Symbol(split(string(key),"_")[1])
        unit = units[var]
        print(io, "    $key: $val $unit \n")
    end
end

# Release property getters
_temperature(r::Release) = r.T
_pressure(r::Release) = r.P
_mass_rate(r::Release) = r.ṁ
_duration(r::Release) = r.Δt
_mass(r::Release) = _mass_rate(r)*_duration(r)
_diameter(r::Release) = r.d
_area(r::Release) = (π/4)*r.d^2
_velocity(r::Release) = r.u
_flowrate(r::Release) = _area(r)*_velocity(r)
_height(r::Release) = r.h
_liquid_fraction(r::Release) = r.f_l

# Stability classes
struct ClassA <: StabilityClass end
struct ClassB <: StabilityClass end
struct ClassC <: StabilityClass end
struct ClassD <: StabilityClass end
struct ClassE <: StabilityClass end
struct ClassF <: StabilityClass end

Base.isapprox(a::Atmosphere, b::Atmosphere) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

function Base.show(io::IO, mime::MIME"text/plain", a::Atmosphere)
    a_type = split(string(typeof(a)),"{")[1]
    print(io, "$a_type atmosphere:\n")
    for key in fieldnames(typeof(a))
        val =  getproperty(a, key)
        unit = units[key]
        print(io, "    $key: $val $unit \n")
    end
end

# Atmosphere property getters
_temperature(a::Atmosphere) = a.T
_pressure(a::Atmosphere) = a.P
_windspeed(a::Atmosphere) = a.u
_windspeed_height(a::Atmosphere) = a.h
_stability(a::Atmosphere) = a.stability
_density(a::Atmosphere) = _density(a, _temperature(a), _pressure(a))

include("simple_atmosphere.jl")

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

# property getters
_lapse_rate(s::Scenario) = _lapse_rate(s.atmosphere)
_atmosphere_temperature(s::Scenario) = _temperature(s.atmosphere)
_release_temperature(s::Scenario) = _temperature(s.release)
_atmosphere_pressure(s::Scenario) = _pressure(s.atmosphere)
_release_pressure(s::Scenario) = _pressure(s.release)
_mass_rate(s::Scenario) = _mass_rate(s.release)
_duration(s::Scenario) = _duration(s.release)
_release_mass(s::Scenario) = _mass(s.release)
_release_diameter(s::Scenario) = _diameter(s.release)
_release_area(s::Scenario) = _area(s.release)
_release_velocity(s::Scenario) = _velocity(s.release)
_release_flowrate(s::Scenario) = _flowrate(s.release)
_release_height(s::Scenario) = _height(s.release)
_release_liquid_fraction(s::Scenario) = _liquid_fraction(s.release)
_windspeed(s::Scenario) = _windspeed(s.atmosphere)
_windspeed(s::Scenario, z::Number, es=DefaultSet) = _windspeed(s.atmosphere, z, es)
_windspeed_height(s::Scenario) = _windspeed_height(s.atmosphere)
_stability(s::Scenario) = _stability(s.atmosphere)
_atmosphere_density(s::Scenario) = _density(s.atmosphere)
_release_density(s::Scenario) = _density(s.substance, _release_liquid_fraction(s), _release_temperature(s), _release_pressure(s))

# Default equation set
struct DefaultSet <: EquationSet end

