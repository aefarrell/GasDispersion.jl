struct Release
    mass_rate::Number
    duration::Number
    diameter::Number
    velocity::Number
    height::Number
    pressure::Number
    temperature::Number
    density::Number
end
Release(; mass_rate, duration, diameter, velocity, height, pressure,
    temperature, density) = Release(mass_rate, duration, diameter,
    velocity, height, pressure, temperature, density)


struct Scenario
    release::Release
    atmosphere::Atmosphere
end
Scenario(;release, atmosphere=Ambient()) = Scenario(release, atmosphere)

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
    :stability => ""
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

Base.isapprox(a::Release, b::Release) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

Base.isapprox(a::Scenario, b::Scenario) = (a.release≈b.release)&&(a.atmosphere≈b.atmosphere)
