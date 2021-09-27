struct Scenario
    mass_emission_rate::Union{Missing,Number}
    release_duration::Union{Missing,Number}
    jet_diameter::Union{Missing,Number}
    jet_velocity::Union{Missing,Number}
    jet_density::Union{Missing,Number}
    release_pressure::Union{Missing,Number}
    release_temperature::Union{Missing,Number}
    release_height::Union{Missing,Number}
    windspeed::Union{Missing,Number}
    ambient_density::Union{Missing,Number}
    ambient_pressure::Union{Missing,Number}
    ambient_temperature::Union{Missing,Number}
    pasquill_gifford::Union{Missing,String}
end

function Scenario(d::AbstractDict)
    return Scenario(
    get(Missing, d, :mass_emission_rate),
    get(Missing, d, :release_duration),
    get(Missing, d, :jet_diameter),
    get(Missing, d, :jet_velocity),
    get(Missing, d, :jet_density),
    get(Missing, d, :release_pressure),
    get(Missing, d, :release_temperature),
    get(Missing, d, :release_height),
    get(Missing, d, :windspeed),
    get(Missing, d, :ambient_density),
    get(Missing, d, :ambient_pressure),
    get(Missing, d, :ambient_temperature),
    get(Missing, d, :pasquill_gifford),
    )
end

function Scenario(base_scenario::Scenario; kwargs...)
    d = Dict([ key => getproperty(base_scenario, key) for key in fieldnames(Scenario)])
    for key in filter(k -> k ∈ fieldnames(Scenario), keys(kwargs))
        d[key] = kwargs[key]
    end
    return Scenario(d)
end

units = Dict{Symbol,String}([
    :mass_emission_rate => "kg/s",
    :release_duration => "s",
    :jet_diameter => "m",
    :jet_velocity => "m/s",
    :jet_density => "kg/m^3",
    :release_pressure => "Pa",
    :release_temperature => "K",
    :release_height => "m",
    :windspeed => "m/s",
    :ambient_density => "kg/m^3",
    :ambient_pressure => "Pa",
    :ambient_temperature => "K",
    :pasquill_gifford => ""
])

function Base.show(io::IO, ::MIME"text/plain", s::Scenario)
    print(io, "Release scenario:\n    ")
    for key in fieldnames(Scenario)
        val =  getproperty(s, key)
        unit = units[key]
        print(io, "$key: $val $unit \n    ")
    end
end
