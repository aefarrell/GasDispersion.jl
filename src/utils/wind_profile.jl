struct DefaultWind <: PowerLawWind end

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
                    F<:Number,S<:StabilityClass,W<:PowerLawWind,SX,SY,SZ}
    u0 = _windspeed(a)
    z0 = _windspeed_height(a)
    return _windspeed(u0,z0,z,S,W)
end

_windspeed(a::Atmosphere,z::Number) = _windspeed(a,z,DefaultSet())


struct BusingerWind <: MoninObukhovWind end
"""
    _windspeed(z, u, zR, L, stability_class, BusingerWind; k=0.35)
returns the windspeed function u(z) for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

# References
+ Businger, J. A., J. C. Wyngaard, Y. Izumi, and E. F. Bradley. 1971. "Flux-Profile Relationships in the Atmospheric Surfaace Layer." *Journal of the Atmospheric Sciences*. 28, 181-189. doi: https://doi.org/10.1175/1520-0469(1971)028<0181:FPRITA>2.0.CO;2 
+ Paulson, C. A. 1970. "The Mathematical Representation of Wind Speed and Temperature Profiles in the Unstable Atmospheric Surface Layer." *Journal of Applied Meteorology and Climatology*. 9(6): 857-861. doi: https://doi.org/10.1175/1520-0450(1970)009<0857:TMROWS>2.0.CO;2

# Arguments
- `u` friction velocity
- `zR` surface roughness
- `L`  Monin-Obukhov length
- `stability_class` Pasquill stability class (A, B, C, D, E, F)
- `k`  von Karman's constant, 0.35
"""
function _windspeed(z::Number, u::Number, zR::Number, L::Number, ::Union{Type{ClassA},Type{ClassB},Type{ClassC}}, ::Type{BusingerWind}; k=0.35)
    x = (1-15*(z/L))^0.25 # Paulson, pg 858; Businger, eqn 8 γ₂=15
    ψ₁ = 2log((1+x)/2) + log((1+x^2)/2) - 2atan(x) + π/2 # Paulson, pg 858
    return (u/k)*(log(z/zR)-ψ₁) # Paulson, equation 2
end

function _windspeed(z::Number, u::Number, zR::Number, L::Number, ::Type{ClassD}, ::Type{BusingerWind}; k=0.35)
    return (u/k)*log(z/zR) # Businger, eqn 22
end

function _windspeed(z::Number, u::Number, zR::Number, L::Number, ::Union{Type{ClassE},Type{ClassF}}, ::Type{BusingerWind}; k=0.35)
    return (u/k)*(log(z/zR) + 4.7*(z/L)) # Businger, eqn 29
end

function _windspeed(a::SimpleAtmosphere{F,S},z::Number,::BasicEquationSet{W,SX,SY,SZ}; k=0.35) where {
                    F<:Number,S<:StabilityClass,W<:MoninObukhovWind,SX,SY,SZ}

    uₐ = _windspeed(a)
    zₐ = _windspeed_height(a)
    zR = _surface_roughness(a)
    L  = _monin_obukhov(zR,S,W)
    u⁺ = k*uₐ/_windspeed(zₐ,1,zR,L,S,W; k=1)
    return _windspeed(z,u⁺,zR,L,S,W; k=k)
end