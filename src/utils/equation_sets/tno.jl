# Correlations from the TNO Yellow Book


struct TNOWind <: MoninObukhovWind end
struct TNOPlumeσy <: DispersionFunction end
struct TNOPlumeσz <: DispersionFunction end
TNOPlume = BasicEquationSet{TNOWind,Nothing,TNOPlumeσy,TNOPlumeσz}

"""
    _monin_obukhov(roughness, StabilityClass, TNOWind)
returns the Monin-Obukhov length for a given Pasquill-Gifford stability class
and surface roughness (in meters)

# References
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.

"""
_monin_obukhov(zR::Number, ::Type{ClassA}, ::TNOWind) = 33.162/log10(min(zR,0.5)/1117)
_monin_obukhov(zR::Number, ::Type{ClassB}, ::TNOWind) = 32.258/log10(min(zR,0.5)/11.46)
_monin_obukhov(zR::Number, ::Type{ClassC}, ::TNOWind) = 51.787/log10(min(zR,0.5)/1.324)
_monin_obukhov(zR::Number, ::Type{ClassD}, ::TNOWind) = Inf
_monin_obukhov(zR::Number, ::Type{ClassE}, ::TNOWind) = -48.330/log10(min(zR,0.5)/1.262)
_monin_obukhov(zR::Number, ::Type{ClassF}, ::TNOWind) = -31.325/log10(min(zR,0.5)/19.36)


"""
    _windspeed(z, u, zR, L, stability_class, TNOWind; k=0.4)
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
function _windspeed(z::Number, u::Number, zR::Number, L::Number, ::Union{Type{ClassA},Type{ClassB},Type{ClassC}}, ::Type{TNOWind}; k=0.4)
    a(z) = (1-16*(z/L))^0.25 # Bakkum calls this Ψ′ but that is awful notation
    Ψ(z) = 2log((1+a(z))/2) + log((1+a(z)^2)/2) - 2atan(a(z)) + π/2
    return (u/k)*(log(z/zR) - Ψ(z) + Ψ(zR)) # Bakkum, eqn 4.32
end

function _windspeed(z::Number, u::Number, zR::Number, L::Number, ::Type{ClassD}, ::Type{TNOWind}; k=0.4)
    return (u/k)*log(z/zR) # Bakkum, eqn 4.32
end

function _windspeed(z::Number, u::Number, zR::Number, L::Number, ::Union{Type{ClassE},Type{ClassF}}, ::Type{TNOWind}; k=0.4)
    return (u/k)*(log(z/zR) + 5*((z-zR)/L)) # Bakkum, eqn 4.32
end

"""
    crosswind_dispersion(x, Plume, StabilityClass, TNO)

Plume crosswind dispersion correlations

# References
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
crosswind_dispersion(x::Number, ::Type{ClassA}, ::Type{TNOPlumeσy}) = 0.527x^0.865
crosswind_dispersion(x::Number, ::Type{ClassB}, ::Type{TNOPlumeσy}) = 0.371x^0.866
crosswind_dispersion(x::Number, ::Type{ClassC}, ::Type{TNOPlumeσy}) = 0.209x^0.897
crosswind_dispersion(x::Number, ::Type{ClassD}, ::Type{TNOPlumeσy}) = 0.128x^0.905
crosswind_dispersion(x::Number, ::Type{ClassE}, ::Type{TNOPlumeσy}) = 0.098x^0.902
crosswind_dispersion(x::Number, ::Type{ClassF}, ::Type{TNOPlumeσy}) = 0.065x^0.902

"""
    vertical_dispersion(x, Plume, StabilityClass, TNO)

Plume vertical dispersion correlations

References:
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
vertical_dispersion(x::Number, ::Type{ClassA}, ::Type{TNOPlumeσz}) = 0.28x^0.90
vertical_dispersion(x::Number, ::Type{ClassB}, ::Type{TNOPlumeσz}) = 0.23x^0.85
vertical_dispersion(x::Number, ::Type{ClassC}, ::Type{TNOPlumeσz}) = 0.22x^0.80
vertical_dispersion(x::Number, ::Type{ClassD}, ::Type{TNOPlumeσz}) = 0.20x^0.76
vertical_dispersion(x::Number, ::Type{ClassE}, ::Type{TNOPlumeσz}) = 0.15x^0.73
vertical_dispersion(x::Number, ::Type{ClassF}, ::Type{TNOPlumeσz}) = 0.12x^0.67

# Puff correlations
struct TNOPuffσx <: DispersionFunction end
struct TNOPuffσy <: DispersionFunction end
struct TNOPuffσz <: DispersionFunction end
TNOPuff = BasicEquationSet{TNOWind,TNOPuffσx,TNOPuffσy,TNOPuffσz}

"""
    crosswind_dispersion(x, Puff, StabilityClass, TNO)

Puff crosswind dispersion correlations

# References
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
crosswind_dispersion(x::Number, stability::Any, ::Type{TNOPuffσy}) = 0.5*crosswind_dispersion(x,stability,TNOPlumeσy)


"""
    vertical_dispersion(x, Puff, StabilityClass, TNO)

Puff vertical dispersion correlations

References:
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
vertical_dispersion(x::Number, ::Type{ClassA}, ::Type{TNOPuffσz}) = 0.28x
vertical_dispersion(x::Number, ::Type{ClassB}, ::Type{TNOPuffσz}) = 0.23x
vertical_dispersion(x::Number, ::Type{ClassC}, ::Type{TNOPuffσz}) = 0.22x
vertical_dispersion(x::Number, ::Type{ClassD}, ::Type{TNOPuffσz}) = 0.20x
vertical_dispersion(x::Number, ::Type{ClassE}, ::Type{TNOPuffσz}) = 0.15x
vertical_dispersion(x::Number, ::Type{ClassF}, ::Type{TNOPuffσz}) = 0.12x

"""
    downwind_dispersion(x, Puff, StabilityClass, TNO)

Puff downwind dispersion correlations

References:
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
downwind_dispersion(x::Number, ::Any, ::Type{TNOPuffσx}) = 0.13x