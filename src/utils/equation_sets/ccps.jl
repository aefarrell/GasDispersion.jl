# Correlations from the CCPS *Guidelines for Consequence Analysis of Chemical Releases*


# Correlations for rural terrain
struct CCPSRural <: EquationSet end

"""
    _windspeed(u0::Number,z0::Number,z::Number, stability, ::CCPSRural)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

Assumes rural terrain.

# References
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)
"""
_windspeed(u0::Number,z0::Number,z::Number,::Union{Type{ClassA},Type{ClassB}},::CCPSRural) = u0*(z/z0)^0.07
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassC},::CCPSRural) = u0*(z/z0)^0.10
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassD},::CCPSRural) = u0*(z/z0)^0.15
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassE},::CCPSRural) = u0*(z/z0)^0.35
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassF},::CCPSRural) = u0*(z/z0)^0.55


"""
    crosswind_dispersion(x, Plume, StabilityClass, ::CCPSRural)

Plume crosswind dispersion correlations, for rural terrain

# References
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)
"""
crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassA}, ::CCPSRural) = 0.22x/√(1+0.0001x)
crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassB}, ::CCPSRural) = 0.16x/√(1+0.0001x)
crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassC}, ::CCPSRural) = 0.11x/√(1+0.0001x)
crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassD}, ::CCPSRural) = 0.08x/√(1+0.0001x)
crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassE}, ::CCPSRural) = 0.06x/√(1+0.0001x)
crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassF}, ::CCPSRural) = 0.04x/√(1+0.0001x)


"""
    vertical_dispersion(x, Plume, StabilityClass, ::CCPSRural)

Plume vertical dispersion correlations, for rural terrain

References:
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)
"""
vertical_dispersion(x, ::Type{Plume}, ::Type{ClassA}, ::CCPSRural) = 0.20x
vertical_dispersion(x, ::Type{Plume}, ::Type{ClassB}, ::CCPSRural) = 0.12x
vertical_dispersion(x, ::Type{Plume}, ::Type{ClassC}, ::CCPSRural) = 0.08x/√(1+0.0002x)
vertical_dispersion(x, ::Type{Plume}, ::Type{ClassD}, ::CCPSRural) = 0.06x/√(1+0.0015x)
vertical_dispersion(x, ::Type{Plume}, ::Type{ClassE}, ::CCPSRural) = 0.03x/(1+0.0003x)
vertical_dispersion(x, ::Type{Plume}, ::Type{ClassF}, ::CCPSRural) = 0.016x/(1+0.0003x)


# Correlations for urban terrain
struct CCPSUrban <: EquationSet end

"""
    _windspeed(u0::Number,z0::Number,z::Number, stability, ::CCPSUrban)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

Assumes urban terrain.

# References
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)
"""
_windspeed(u0::Number,z0::Number,z::Number,::Union{Type{ClassA},Type{ClassB}},::CCPSUrban) = u0*(z/z0)^0.15
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassC},::CCPSUrban) = u0*(z/z0)^0.20
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassD},::CCPSUrban) = u0*(z/z0)^0.25
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassE},::CCPSUrban) = u0*(z/z0)^0.40
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassF},::CCPSUrban) = u0*(z/z0)^0.60


"""
    crosswind_dispersion(x, Plume, StabilityClass, ::CCPSUrban)

Plume crosswind dispersion correlations, for urban terrain

# References
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)
"""
crosswind_dispersion(x, ::Type{Plume}, ::Union{Type{ClassA},Type{ClassB}}, ::CCPSUrban) = 0.32x/√(1+0.0004x)
crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassC}, ::CCPSUrban) = 0.22x/√(1+0.0004x)
crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassD}, ::CCPSUrban) = 0.16x/√(1+0.0004x)
crosswind_dispersion(x, ::Type{Plume}, ::Union{Type{ClassE},Type{ClassF}}, ::CCPSUrban) = 0.11x/√(1+0.0004x)


"""
    vertical_dispersion(x, Plume, StabilityClass, ::CCPSUrban)

Plume vertical dispersion correlations, for urban terrain

References:
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)
"""
vertical_dispersion(x, ::Type{Plume}, ::Union{Type{ClassA},Type{ClassB}}, ::CCPSUrban) = 0.24x/√(1+0.001x)
vertical_dispersion(x, ::Type{Plume}, ::Type{ClassC}, ::CCPSUrban) = 0.20x
vertical_dispersion(x, ::Type{Plume}, ::Type{ClassD}, ::CCPSUrban) = 0.14x/√(1+0.003x)
vertical_dispersion(x, ::Type{Plume}, ::Union{Type{ClassE},Type{ClassF}}, ::CCPSUrban) = 0.08x/√(1+0.0015x)


# Correlations for both

"""
    _windspeed(a::Atmosphere, z::Number, ::Union{CCPSUrban, CCPSRural})
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class using a power law correlation.

u = u₀ (z/z₀)ᵖ

Where *p* is a function of the atmospheric stability class and the terrain (urban or rural)

# References
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)
"""
function _windspeed(a::Atmosphere,z::Number,es::Union{CCPSUrban, CCPSRural})
    stab = _stability(a)
    u0 = _windspeed(a)
    z0 = _windspeed_height(a)
    return _windspeed(u0,z0,z,stab,es)
end

# Puff dispersion correlations
"""
    crosswind_dispersion(x, Puff, StabilityClass. ::Union{CCPSUrban,CCPSRural})

Puff crosswind dispersion correlations

References:
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)
"""
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassA}, ::Union{CCPSUrban,CCPSRural}) = 0.18*x^0.92
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassB}, ::Union{CCPSUrban,CCPSRural}) = 0.14*x^0.92
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassC}, ::Union{CCPSUrban,CCPSRural}) = 0.10*x^0.92
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassD}, ::Union{CCPSUrban,CCPSRural}) = 0.06*x^0.92
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassE}, ::Union{CCPSUrban,CCPSRural}) = 0.04*x^0.92
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassF}, ::Union{CCPSUrban,CCPSRural}) = 0.02*x^0.89


"""
    downwind_dispersion(x, Puff, StabilityClass)

Puff downwind dispersion correlations

References:
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)
"""
downwind_dispersion(x, ::Type{Puff}, 
    stab::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD},Type{ClassE},Type{ClassF}}, 
    es::Union{CCPSUrban,CCPSRural}) =crosswind_dispersion(x, Puff, stab, es)


"""
    vertical_dispersion(x, Puff, StabilityClass)

Puff vertical dispersion correlations

References:
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)
"""
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassA}, ::Union{CCPSUrban,CCPSRural}) = 0.60*x^0.75
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassB}, ::Union{CCPSUrban,CCPSRural}) = 0.53*x^0.73
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassC}, ::Union{CCPSUrban,CCPSRural}) = 0.34*x^0.71
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassD}, ::Union{CCPSUrban,CCPSRural}) = 0.15*x^0.70
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassE}, ::Union{CCPSUrban,CCPSRural}) = 0.10*x^0.65
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassF}, ::Union{CCPSUrban,CCPSRural}) = 0.05*x^0.61
