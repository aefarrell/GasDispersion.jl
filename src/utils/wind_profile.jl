struct DefaultWind <: PowerLawWind end

# Power law correlations
windspeed(u0::Number,z0::Number,z::Number,::ClassA,::DefaultWind) = u0*(z/z0)^0.108
windspeed(u0::Number,z0::Number,z::Number,::ClassB,::DefaultWind) = u0*(z/z0)^0.112
windspeed(u0::Number,z0::Number,z::Number,::ClassC,::DefaultWind) = u0*(z/z0)^0.120
windspeed(u0::Number,z0::Number,z::Number,::ClassD,::DefaultWind) = u0*(z/z0)^0.142
windspeed(u0::Number,z0::Number,z::Number,::ClassE,::DefaultWind) = u0*(z/z0)^0.203
windspeed(u0::Number,z0::Number,z::Number,::ClassF,::DefaultWind) = u0*(z/z0)^0.253


"""
    windspeed(a::Atmosphere, z::Number)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

# References

"""
function windspeed(a::Atmosphere,z::Number,es::BasicEquationSet{W,SX,SY,SZ}) where {W<:PowerLawWind,SX,SY,SZ}
    u0 = windspeed(a)
    z0 = _windspeed_height(a)
    stab = _stability(a)
    wind = _wind_equation(es)
    return windspeed(u0,z0,z,stab,wind)
end

windspeed(a::Atmosphere,z::Number) = windspeed(a,z,DefaultSet)

# Monin-Obukhov similarity theory wind profiles
struct BusingerWind <: MoninObukhovWind end
"""
    windspeed(z, u, zR, L, stability_class, BusingerWind; k=0.35)
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
function windspeed(z::Number, u::Number, zR::Number, L::Number, ::UnstableClass, ::BusingerWind; k=0.35)
    x = (1-15*(z/L))^0.25 # Paulson, pg 858; Businger, eqn 8 γ₂=15
    ψ₁ = 2log((1+x)/2) + log((1+x^2)/2) - 2atan(x) + π/2 # Paulson, pg 858
    return (u/k)*(log(z/zR)-ψ₁) # Paulson, equation 2
end

function windspeed(z::Number, u::Number, zR::Number, L::Number, ::NeutralClass, ::BusingerWind; k=0.35)
    return (u/k)*log(z/zR) # Businger, eqn 22
end

function windspeed(z::Number, u::Number, zR::Number, L::Number, ::StableClass, ::BusingerWind; k=0.35)
    return (u/k)*(log(z/zR) + 4.7*(z/L)) # Businger, eqn 29
end

function windspeed(a::Atmosphere,z::Number,es::BasicEquationSet{W,SX,SY,SZ}; k=0.35) where {W<:MoninObukhovWind,SX,SY,SZ}

    stab  = _stability(a)
    wind = _wind_equation(es)
    zR = _surface_roughness(a)
    L  = monin_obuknov(a,es)
    u⁺ = friction_velocity(a,es; k=k)
    return windspeed(z,u⁺,zR,L,stab,wind; k=k)
end

# friction velocity
friction_velocity(a::Atmosphere) = friction_velocity(a,DefaultSet)

function friction_velocity(a::SimpleAtmosphere,es::BasicEquationSet{W,SX,SY,SZ}; k=0.35) where {W<:PowerLawWind,SX,SY,SZ}
    return 0.1*windspeed(a,10,es)
end

function friction_velocity(a::SimpleAtmosphere,es::BasicEquationSet{W,SX,SY,SZ}; k=0.35) where {W<:MoninObukhovWind,SX,SY,SZ}
    uₐ = windspeed(a)
    zₐ = _windspeed_height(a)
    zR = _surface_roughness(a)
    stab = _stability(a)
    wind = _wind_equation(es)
    L  = monin_obuknov(zR,stab,wind)
    return k*uₐ/windspeed(zₐ,1,zR,L,stab,wind; k=1)
end