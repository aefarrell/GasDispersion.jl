# Power law correlations
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassA},::Type{DefaultWind}) = u0*(z/z0)^0.108
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassB},::Type{DefaultWind}) = u0*(z/z0)^0.112
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassC},::Type{DefaultWind}) = u0*(z/z0)^0.120
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassD},::Type{DefaultWind}) = u0*(z/z0)^0.142
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassE},::Type{DefaultWind}) = u0*(z/z0)^0.203
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassF},::Type{DefaultWind}) = u0*(z/z0)^0.253


"""
    _windspeed(a::Atmosphere, z::Number)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

# References

"""
function _windspeed(a::SimpleAtmosphere{F,S},z::Number,::BasicEquationSet{W,SX,SY,SZ}) where {
                    F<:Number,S<:StabilityClass,W<:PowerLawWind,SX<:Any,SY<:Any,SZ<:Any}
    u0 = _windspeed(a)
    z0 = _windspeed_height(a)
    return _windspeed(u0,z0,z,S,W)
end

_windspeed(a::Atmosphere,z::Number) = _windspeed(a,z,DefaultSet())


"""
    _windspeed(z, u, zR, λ, stability_class; k=0.35)
returns the windspeed function u(z) for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

# References
+ Businger et al. 1971

# Arguments
- `u` friction velocity
- `zR` surface roughness
- `λ`  Monin-Obukhov length
- `stability_class` Pasquill stability class (A, B, C, D, E, F)
- `k`  von Karman's constant, 0.35
"""
function _windspeed(z::Number, u::Number, zR::Number, λ::Number, ::Union{Type{ClassA},Type{ClassB},Type{ClassC}}; k=0.35)
    a = (1-15*(z/λ))^0.25
    Ψ = 2*log((1+a)/2) + log((1+a^2)/2) - 2*atan(a) + π/2
    return (u/k)*(log((z+zR)/zR) - Ψ)
end

function _windspeed(z::Number, u::Number, zR::Number, λ::Number, ::Type{ClassD}; k=0.35)
    return (u/k)*log((z+zR)/zR)
end

function _windspeed(z::Number, u::Number, zR::Number, λ::Number, ::Union{Type{ClassE},Type{ClassF}}; k=0.35)
    return (u/k)*(log((z+zR)/zR) - 4.7*(z/λ))
end
