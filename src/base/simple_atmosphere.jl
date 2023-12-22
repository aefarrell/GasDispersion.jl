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

_lapse_rate(a::SimpleAtmosphere{<:Number,ClassE}) = 0.020
_lapse_rate(a::SimpleAtmosphere{<:Number,ClassF}) = 0.035
_lapse_rate(a::SimpleAtmosphere{<:Number,<:StabilityClass}) = nothing
_density(a::SimpleAtmosphere, T, P) = P/(a.Rs*T)
_rel_humidity(a::SimpleAtmosphere) = a.rh