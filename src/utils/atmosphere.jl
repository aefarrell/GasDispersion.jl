struct Ambient <: Atmosphere
    pressure::Number
    temperature::Number
    density::Number
    windspeed::Number
    stability::String
end
Ambient(; pressure=101325, temperature=298.15, density=1.225, windspeed=1.5,
          stability="F") = Ambient(pressure,temperature,density,windspeed,stability)
