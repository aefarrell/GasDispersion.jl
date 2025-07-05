# Correlations from the TNO Yellow Book


struct TNOWind <: MoninObukhovWind end
struct TNOPlumeσy <: DispersionFunction end
struct TNOPlumeσz <: DispersionFunction end
const TNOPlume = BasicEquationSet(TNOWind(),nothing,TNOPlumeσy(),TNOPlumeσz())

# Default surface roughness for TNO models is 0.1
function windspeed(a::SimpleAtmosphere,z::Number,es::BasicEquationSet{W,SX,SY,SZ}; k=0.4) where {W<:TNOWind,SX,SY,SZ}
    zR = 0.1
    stab = _stability(a)
    L  = monin_obuknov(zR, stab, TNOWind())
    uₐ = windspeed(a)
    zₐ = _windspeed_height(a)
    u⁺ = k*uₐ/windspeed(zₐ, 1, zR, L, stab, TNOWind(); k=1)
    return windspeed(z, u⁺, zR, L, stab, TNOWind(); k=k)
end

"""
    monin_obuknov(roughness, StabilityClass, TNOWind)
returns the Monin-Obukhov length for a given Pasquill-Gifford stability class
and surface roughness (in meters)

# References
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.

"""
monin_obuknov(zR::Number, ::ClassA, ::TNOWind) = 33.162/log10(min(zR,0.5)/1117)
monin_obuknov(zR::Number, ::ClassB, ::TNOWind) = 32.258/log10(min(zR,0.5)/11.46)
monin_obuknov(zR::Number, ::ClassC, ::TNOWind) = 51.787/log10(min(zR,0.5)/1.324)
monin_obuknov(zR::Number, ::ClassD, ::TNOWind) = Inf
monin_obuknov(zR::Number, ::ClassE, ::TNOWind) = -48.330/log10(min(zR,0.5)/1.262)
monin_obuknov(zR::Number, ::ClassF, ::TNOWind) = -31.325/log10(min(zR,0.5)/19.36)


"""
    windspeed(z, u, zR, L, stability_class, TNOWind; k=0.4)
returns the windspeed function u(z) for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

# References
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.

# Arguments
- `u` friction velocity
- `zR` surface roughness
- `L`  Monin-Obukhov length
- `stability_class` Pasquill stability class (A, B, C, D, E, F)
- `k`  von Karman's constant, 0.40
"""
function windspeed(z::Number, u::Number, zR::Number, L::Number, ::UnstableClass, ::TNOWind; k=0.4)
    a(z) = (1-16*(z/L))^0.25 # Bakkum calls this Ψ′ but that is awful notation
    Ψ(z) = 2log((1+a(z))/2) + log((1+a(z)^2)/2) - 2atan(a(z)) + π/2
    return (u/k)*(log(z/zR) - Ψ(z) + Ψ(zR)) # Bakkum, eqn 4.32
end

function windspeed(z::Number, u::Number, zR::Number, L::Number, ::NeutralClass, ::TNOWind; k=0.4)
    return (u/k)*log(z/zR) # Bakkum, eqn 4.32
end

function windspeed(z::Number, u::Number, zR::Number, L::Number, ::StableClass, ::TNOWind; k=0.4)
    return (u/k)*(log(z/zR) + 5*((z-zR)/L)) # Bakkum, eqn 4.32
end

"""
    crosswind_dispersion(x, StabilityClass, TNOPlumeσy)

Plume crosswind dispersion correlations

# References
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
crosswind_dispersion(x::Number, ::ClassA, ::TNOPlumeσy) = 0.527x^0.865
crosswind_dispersion(x::Number, ::ClassB, ::TNOPlumeσy) = 0.371x^0.866
crosswind_dispersion(x::Number, ::ClassC, ::TNOPlumeσy) = 0.209x^0.897
crosswind_dispersion(x::Number, ::ClassD, ::TNOPlumeσy) = 0.128x^0.905
crosswind_dispersion(x::Number, ::ClassE, ::TNOPlumeσy) = 0.098x^0.902
crosswind_dispersion(x::Number, ::ClassF, ::TNOPlumeσy) = 0.065x^0.902

"""
    vertical_dispersion(x, StabilityClass, TNOPlumeσz)

Plume vertical dispersion correlations

References:
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
vertical_dispersion(x::Number, ::ClassA, ::TNOPlumeσz) = 0.28x^0.90
vertical_dispersion(x::Number, ::ClassB, ::TNOPlumeσz) = 0.23x^0.85
vertical_dispersion(x::Number, ::ClassC, ::TNOPlumeσz) = 0.22x^0.80
vertical_dispersion(x::Number, ::ClassD, ::TNOPlumeσz) = 0.20x^0.76
vertical_dispersion(x::Number, ::ClassE, ::TNOPlumeσz) = 0.15x^0.73
vertical_dispersion(x::Number, ::ClassF, ::TNOPlumeσz) = 0.12x^0.67

# Puff correlations
struct TNOPuffσx <: DispersionFunction end
struct TNOPuffσy <: DispersionFunction end
struct TNOPuffσz <: DispersionFunction end
const TNOPuff = BasicEquationSet(TNOWind(),TNOPuffσx(),TNOPuffσy(),TNOPuffσz())

"""
    crosswind_dispersion(x, StabilityClass, TNOPuffσy)

Puff crosswind dispersion correlations

# References
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
crosswind_dispersion(x::Number, stability::StabilityClass, ::TNOPuffσy) = 0.5*crosswind_dispersion(x,stability,TNOPlumeσy())


"""
    vertical_dispersion(x, StabilityClass, TNOPuffσz)

Puff vertical dispersion correlations

References:
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
vertical_dispersion(x::Number, ::ClassA, ::TNOPuffσz) = 0.28x
vertical_dispersion(x::Number, ::ClassB, ::TNOPuffσz) = 0.23x
vertical_dispersion(x::Number, ::ClassC, ::TNOPuffσz) = 0.22x
vertical_dispersion(x::Number, ::ClassD, ::TNOPuffσz) = 0.20x
vertical_dispersion(x::Number, ::ClassE, ::TNOPuffσz) = 0.15x
vertical_dispersion(x::Number, ::ClassF, ::TNOPuffσz) = 0.12x

"""
    downwind_dispersion(x, StabilityClass, TNOPuffσx)

Puff downwind dispersion correlations

References:
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
downwind_dispersion(x::Number, ::StabilityClass, ::TNOPuffσx) = 0.13x