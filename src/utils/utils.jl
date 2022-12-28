# methods to get properties
_mass_rate(s::Scenario) = s.release.mass_rate
_duration(s::Scenario) = s.release.duration
_release_diameter(s::Scenario) = s.release.diameter
_release_area(s::Scenario) = (π/4)*s.release.diameter^2
_release_velocity(s::Scenario) = s.release.velocity
_release_flowrate(s::Scenario) = _release_area(s)*_release_velocity(s)
_release_mass(s::Scenario) = s.release.mass_rate*s.release.duration
_release_height(s::Scenario) = s.release.height
_release_pressure(s::Scenario) = s.release.pressure
_release_temperature(s::Scenario) = s.release.temperature
_release_density(s::Scenario) = _release_density(s.release.fraction_liquid,
                                                 s.release.temperature,
                                                 s.substance)
_atmosphere_pressure(s::Scenario) = s.atmosphere.pressure
_atmosphere_temperature(s::Scenario) = s.atmosphere.temperature
_atmosphere_density(s::Scenario) = s.atmosphere.density
_windspeed(s::Scenario) = s.atmosphere.windspeed
_stability(s::Scenario) = s.atmosphere.stability

function _release_density(f_l, T, s::Substance{<:Number,<:Union{Function, Number}})
    f_g = 1 - f_l
    ρ_l = s.ρ_l
    ρ_g = s.ρ_g

    return 1/(f_l/ρ_l + f_g/ρ_g)
end

function _release_density(f_l, T, s::Substance{<:Function,<:Union{Function, Number}})
    f_g = 1 - f_l
    ρ_l = s.ρ_l(T)
    ρ_g = s.ρ_g(T)

    return 1/(f_l/ρ_l + f_g/ρ_g)
end

# method to print the scenario to the screen nicely
units = Dict{Symbol,String}([
    :mass_rate => "kg/s",
    :duration => "s",
    :diameter => "m",
    :velocity => "m/s",
    :height => "m",
    :density => "kg/m^3",
    :pressure => "Pa",
    :temperature => "K",
    :windspeed => "m/s",
    :stability => "",
    :fraction_liquid => ""
])

function Base.show(io::IO, ::MIME"text/plain", s::Scenario)
    r = s.release
    a = s.atmosphere
    print(io, "Atmospheric conditions:\n    ")
    for key in fieldnames(typeof(a))
        val =  getproperty(a, key)
        unit = units[key]
        print(io, "$key: $val $unit \n    ")
    end
    print(io, "\nRelease conditions:\n    ")
    for key in fieldnames(typeof(r))
        val =  getproperty(r, key)
        unit = units[key]
        print(io, "$key: $val $unit \n    ")
    end
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
