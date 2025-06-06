# units for the show function
const UNITS = Dict{Symbol,String}([
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
                                    :k => "",
                                    :name => "",
                                    :MW => "kg/mol",
                                    :Δh => "J/kg",
                                    :Cp => "J/kg/K"
                                ])

include("substance_correls.jl")

const CALLABLE = Union{Function,PropertyCorrelation}
const CC = Union{Number,CALLABLE}

"""
    Substance{N,VAP,D_G,D_L,F,H,CP_G,CP_L}(;kwargs...)

A simple container for the physical and thermal properties of substances.

`vapor_pressure`, `latent_heat`, `gas_heat_capacity`, and `liquid_heat_capacity` can be functions of temperature (in Kelvin) or constants.

`gas_density` and `liquid_density` can be functions of temperature (in Kelvin) and pressure (in Pascal) or constants.

# Arguments
- `name<:Union{AbstractString,Symbol}`: the name of the substance
- `molar_weight::Number`: the molar weight, kg/mol
- `vapor_pressure<:Union{Number,Function,Nothing}=nothing`: the vapor pressure, Pa. If `nothing` then the Clausius-Clapeyron equation is used to derive a vapor pressure curve
- `gas_density<:Union{Number,Function,Nothing}=nothing`: the gas density, kg/m³. If `nothing` the ideal gas law is used.
- `liquid_density<:Union{Number,Function}`: the liquid density, kg/m³.
- `reference_temp::Number=288.15`: the reference temperature for the given properties, K.
- `reference_pressure::Number=101325`: the reference pressure for the given properties, Pa.
- `k::Number=1.4`: the isentropic expansion factor, cp/cv, unitless.
- `boiling_temp<:Union{Number,Function}::Number`: the normal boiling temperature, K.
- `latent_heat<:Union{Number,Function}`: the latent heat of vaporization, J/kg.
- `gas_heat_capacity<:Union{Number,Function}`: the gas heat capacity, J/kg/K.
- `liquid_heat_capacity<:Union{Number,Function}`: the liquid heat capacity, J/kg/K.

"""
struct Substance{N<:Union{AbstractString,Symbol},VAP<:CC,D_G<:CC,D_L<:CC,F<:Number,H<:CC,CP_G<:CC,CP_L<:CC} <: AbstractSubstance
    name::N
    MW::F  # molar weight, kg/mol
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
   
    # initialize with ideal gas
    if isnothing(gas_density)
        gas_density = reference_pressure*molar_weight/(R_GAS_CONST*reference_temp)
    end

    # if no vapour pressure correlation given, use Clapeyron equation
    if isnothing(vapor_pressure)
        if latent_heat isa Number
            Δh = latent_heat
        else
            Δh = latent_heat(boiling_temp)
        end

        B = Δh*molar_weight/R_GAS_CONST
        A = B/(boiling_temp)
        C = 0.0
        vapor_pressure = Antoine(A,B,C)
    end

    molar_weight,reference_temp,reference_pressure,k,boiling_temp = promote(molar_weight,reference_temp,reference_pressure,k,boiling_temp)

    return Substance(name,molar_weight,vapor_pressure,gas_density,liquid_density,
    reference_temp,reference_pressure,k,boiling_temp,latent_heat,
    gas_heat_capacity,liquid_heat_capacity)
end

Substance(;name,molar_weight,vapor_pressure=nothing,gas_density=nothing,liquid_density,reference_temp=288.15,
          reference_pressure=101325,k=1.4,boiling_temp,latent_heat,gas_heat_capacity,
          liquid_heat_capacity) = Substance(name,molar_weight,vapor_pressure,gas_density,liquid_density,
          reference_temp,reference_pressure,k,boiling_temp,latent_heat,gas_heat_capacity,
          liquid_heat_capacity)

Base.isapprox(a::Substance, b::Substance) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

function Base.show(io::IO, mime::MIME"text/plain", s::S) where { S<:Substance}
    s_type = split(string(S),"{")[1]
    print(io, "$s_type: $(s.name) \n")
    for key in fieldnames(S) 
        if key != :name
            val =  getproperty(s, key)
            var =  Symbol(split(string(key),"_")[1])
            unit = UNITS[var]
            print(io, "    $key: $val $unit \n")
        end
    end
end

# Substance property getters
_MW(s::Substance) = s.MW
_vapor_pressure(s::Substance, T) = s.P_v(T)

_cp_gas(s::Substance{<:Any,<:Any,<:Any,<:Any,<:Any,<:Number,<:Any},T) = s.Cp_g
_cp_gas(s::Substance{<:Any,<:Any,<:Any,<:Any,<:Any,<:CALLABLE,<:Any},T) = s.Cp_g(T)
_cp_gas(s::Substance) = _cp_gas(s, s.T_ref)

_cp_liquid(s::Substance{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Number},T) = s.Cp_l
_cp_liquid(s::Substance{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:CALLABLE},T) = s.Cp_l(T)
_cp_liquid(s::Substance) = _cp_liquid(s, s.T_ref)

_boiling_temperature(s::Substance) = s.T_b
_latent_heat(s::Substance{<:Any,<:Any,<:Any,<:Any,<:Any,<:Number},T) = s.Δh_v
_latent_heat(s::Substance{<:Any,<:Any,<:Any,<:Any,<:Any,<:CALLABLE},T) = s.Δh_v(T)
_latent_heat(s::Substance) = _latent_heat(s, s.T_b)

# density functions
_liquid_density(s::Substance{<:Any,<:Any,<:Any,<:Number}, T::Number, P::Number) = s.ρ_l
_liquid_density(s::Substance{<:Any,<:Any,<:Any,<:CALLABLE}, T::Number, P::Number) = s.ρ_l(T,P)
_liquid_density(s::Substance) = _liquid_density(s, s.T_ref, s.P_ref)

_gas_density(s::Substance{<:Any,<:Any,<:Number}, T::Number, P::Number) = s.ρ_g*(s.T_ref/T)*(P/s.P_ref)
_gas_density(s::Substance{<:Any,<:Any,<:CALLABLE}, T::Number, P::Number) = s.ρ_g(T,P)
_gas_density(s::Substance) = _gas_density(s, s.T_ref, s.P_ref)

function _density(s::Substance, f_l, T, P)
    f_g = 1 - f_l
    ρ_l = _liquid_density(s,T,P)
    ρ_g = _gas_density(s,T,P)

    return 1/(f_l/ρ_l + f_g/ρ_g)
end



"""
    HorizontalJet{<:Number}(;kwargs...)<:Release

A simple container for parameters of a horizontal jet release.

# Arguments
- `mass_rate::Number`: the mass emission rate of the substance, kg/s.
- `duration::Number`: the duration of the release, s.
- `diameter::Number`: the diameter of the jet, m.
- `velocity::Number`: the average velocity of the jet, m/s.
- `height::Number`: the height of the jet center, m.
- `pressure::Number`: the pressure at the jet exit, Pa.
- `temperature::Number`: the temperature at the jet exit, K.
- `fraction_liquid::Number`: the fraction of the release that is liquid (for gas liquid mixtures), vol fraction.

"""
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

"""
    VerticalJet{<:Number}(;kwargs...)<:Release

A simple container for parameters of a vertical jet release.

# Arguments
- `mass_rate::Number`: the mass emission rate of the substance, kg/s.
- `duration::Number`: the duration of the release, s.
- `diameter::Number`: the diameter of the jet, m.
- `velocity::Number`: the average velocity of the jet, m/s.
- `height::Number`: the height of the jet center, m.
- `pressure::Number`: the pressure at the jet exit, Pa.
- `temperature::Number`: the temperature at the jet exit, K.
- `fraction_liquid::Number`: the fraction of the release that is liquid (for gas liquid mixtures), vol fraction.

"""
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
        unit = UNITS[var]
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
        unit = UNITS[key]
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

"""
    Scenario{<:Substance,<:Release,<:Atmosphere}(s<:Substance,r<:Release,a<:Atmosphere=SimpleAtmosphere())

A chemical release scenario.

"""
struct Scenario{S<:AbstractSubstance,R<:Release,A<:Atmosphere}
    substance::S
    release::R
    atmosphere::A
end
Scenario(substance,release) = Scenario(substance,release,SimpleAtmosphere())
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
_windspeed(s::Scenario, z::Number, es=DefaultSet()) = _windspeed(s.atmosphere, z, es)
_windspeed_height(s::Scenario) = _windspeed_height(s.atmosphere)
_stability(s::Scenario) = _stability(s.atmosphere)
_atmosphere_density(s::Scenario) = _density(s.atmosphere)
_release_density(s::Scenario) = _density(s.substance, _release_liquid_fraction(s), _release_temperature(s), _release_pressure(s))

# Default equation set
struct BasicEquationSet{WIND,SIGMAX,SIGMAY,SIGMAZ} <: EquationSet end