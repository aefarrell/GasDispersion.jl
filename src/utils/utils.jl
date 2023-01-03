# methods to get properties
_temperature(a::Atmosphere) = a.temperature
_temperature(r::Release) = r.temperature
_atmosphere_temperature(s::Scenario) = _temperature(s.atmosphere)
_release_temperature(s::Scenario) = _temperature(s.release)

_pressure(a::Atmosphere) = a.pressure
_pressure(r::Release) = r.pressure
_atmosphere_pressure(s::Scenario) = _pressure(s.atmosphere)
_release_pressure(s::Scenario) = _pressure(s.release)

_mass_rate(s::Scenario) = s.release.mass_rate
_duration(s::Scenario) = s.release.duration
_release_diameter(s::Scenario) = s.release.diameter
_release_area(s::Scenario) = (π/4)*_release_diameter(s)^2
_release_velocity(s::Scenario) = s.release.velocity
_release_flowrate(s::Scenario) = _release_area(s)*_release_velocity(s)
_release_mass(s::Scenario) = _mass_rate(s)*_duration(s)
_release_height(s::Scenario) = s.release.height
_release_liquid_fraction(s::Scenario) = s.release.fraction_liquid

_windspeed(a::Atmosphere) = a.windspeed
_windspeed(s::Scenario) = _windspeed(s.atmosphere)
_windspeed(s::Scenario, z) = _windspeed(s.atmosphere, z)
_windspeed_height(a::Atmosphere) = a.windspeed_height
_windspeed_height(s::Scenario) = _windspeed_height(s.atmosphere)
_stability(a::Atmosphere) = a.stability
_stability(s::Scenario) = _stability(s.atmosphere)

# density functions
_liquid_density(s::Substance) = _liquid_density(s, s.T_ref, s.P_ref)
_liquid_density(s::Substance{<:Number,<:Union{Function, Number}}, T::Number, P::Number) = s.ρ_l
_liquid_density(s::Substance{<:Function,<:Union{Function, Number}}, T::Number, P::Number) = s.ρ_l(T,P)

_gas_density(s::Substance) = _gas_density(s, s.T_ref, s.P_ref)
_gas_density(s::Substance{<:Number,<:Union{Function, Number}}, T::Number, P::Number) = s.ρ_g*(s.T_ref/T)*(P/s.P_ref)
_gas_density(s::Substance{<:Function,<:Union{Function, Number}}, T::Number, P::Number) = s.ρ_g(T,P)

function _density(s::Substance, f_l, T, P)
    f_g = 1 - f_l
    ρ_l = _liquid_density(s,T,P)
    ρ_g = _gas_density(s,T,P)

    return 1/(f_l/ρ_l + f_g/ρ_g)
end

_density(a::Atmosphere, T, P) = a.density
_density(a::Atmosphere) = _density(a, _temperature(a), _pressure(a))
_atmosphere_density(s::Scenario) = _density(s.atmosphere)
_release_density(s::Scenario) = _density(s.substance, _release_liquid_fraction(s), _release_temperature(s), _release_pressure(s))


# method to print the scenario to the screen nicely
units = Dict{Symbol,String}([
    :mass_rate => "kg/s",
    :duration => "s",
    :diameter => "m",
    :velocity => "m/s",
    :height => "m",
    :windspeed_height => "m",
    :density => "kg/m^3",
    :pressure => "Pa",
    :temperature => "K",
    :windspeed => "m/s",
    :stability => "",
    :fraction_liquid => ""
])

function Base.show(io::IO, mime::MIME"text/plain", s::Substance)
    print(io, "Substance: $(s.name) \n")
end

function Base.show(io::IO, mime::MIME"text/plain", r::Release)
    print(io, "Release conditions:\n")
    for key in fieldnames(typeof(r))
        val =  getproperty(r, key)
        unit = units[key]
        print(io, "    $key: $val $unit \n")
    end
end

function Base.show(io::IO, mime::MIME"text/plain", a::Atmosphere)
    print(io, "Atmospheric conditions:\n")
    for key in fieldnames(typeof(a))
        val =  getproperty(a, key)
        unit = units[key]
        print(io, "    $key: $val $unit \n")
    end
end

function Base.show(io::IO, mime::MIME"text/plain", s::Scenario)
    show(io,mime,s.substance)
    show(io,mime,s.release)
    show(io,mime,s.atmosphere)
end

Base.isapprox(a::Atmosphere, b::Atmosphere) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

Base.isapprox(a::Substance, b::Substance) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

Base.isapprox(a::Release, b::Release) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

Base.isapprox(a::Scenario, b::Scenario) = all( [ a.substance≈b.substance,
                                                 a.release≈b.release,
                                                 a.atmosphere≈b.atmosphere ])


# Pasquill-Gifford dispersion correlations
include("pasquill_gifford.jl")

# correlations for Monin-Obukhov length
include("monin_obukhov.jl")

# correlations for windspeed
include("wind_profile.jl")
