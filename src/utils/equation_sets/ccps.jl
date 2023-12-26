# Correlations from the CCPS *Guidelines for Consequence Analysis of Chemical Releases*


# Correlations for rural terrain
struct CCPSRural <: EquationSet end

"""
    _windspeed(u0::Number,z0::Number,z::Number, stability, CCPSRural)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

Assumes rural terrain.

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
"""
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassA},::Type{CCPSRural}) = u0*(z/z0)^0.07
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassB},::Type{CCPSRural}) = u0*(z/z0)^0.07
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassC},::Type{CCPSRural}) = u0*(z/z0)^0.10
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassD},::Type{CCPSRural}) = u0*(z/z0)^0.15
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassE},::Type{CCPSRural}) = u0*(z/z0)^0.35
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassF},::Type{CCPSRural}) = u0*(z/z0)^0.55


"""
    crosswind_dispersion(x, Plume, StabilityClass, CCPSRural)

Plume crosswind dispersion correlations, for rural terrain

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
+ Briggs, Gary A. 1973. *Diffusion Estimation for Small Emissions. Preliminary Report*. United States. https://doi.org/10.2172/5118833
+ Griffiths, R. F. 1994. "Errors in the use of the Briggs parameterization for atmospheric dispersion coefficients." *Atmospheric Environment* 28(17):2861-2865 https://doi.org/10.1016/1352-2310(94)90086-8
"""
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::Type{CCPSRural}) = 0.22x/√(1+0.0001x)
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::Type{CCPSRural}) = 0.16x/√(1+0.0001x)
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::Type{CCPSRural}) = 0.11x/√(1+0.0001x)
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::Type{CCPSRural}) = 0.08x/√(1+0.0001x)
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::Type{CCPSRural}) = 0.06x/√(1+0.0001x)
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::Type{CCPSRural}) = 0.04x/√(1+0.0001x)


"""
    vertical_dispersion(x, Plume, StabilityClass, CCPSRural)

Plume vertical dispersion correlations, for rural terrain

References:
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
+ Briggs, Gary A. 1973. *Diffusion Estimation for Small Emissions. Preliminary Report*. United States. https://doi.org/10.2172/5118833
+ Griffiths, R. F. 1994. "Errors in the use of the Briggs parameterization for atmospheric dispersion coefficients." *Atmospheric Environment* 28(17):2861-2865 https://doi.org/10.1016/1352-2310(94)90086-8
"""
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::Type{CCPSRural}) = 0.20x
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::Type{CCPSRural}) = 0.12x
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::Type{CCPSRural}) = 0.08x/√(1+0.0002x)
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::Type{CCPSRural}) = 0.06x/√(1+0.0015x)
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::Type{CCPSRural}) = 0.03x/(1+0.0003x)
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::Type{CCPSRural}) = 0.016x/(1+0.0003x)


# Correlations for urban terrain
struct CCPSUrban <: EquationSet end

"""
    _windspeed(u0::Number,z0::Number,z::Number, stability, CCPSUrban)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

Assumes urban terrain.

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
"""
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassA},::Type{CCPSUrban}) = u0*(z/z0)^0.15
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassB},::Type{CCPSUrban}) = u0*(z/z0)^0.15
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassC},::Type{CCPSUrban}) = u0*(z/z0)^0.20
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassD},::Type{CCPSUrban}) = u0*(z/z0)^0.25
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassE},::Type{CCPSUrban}) = u0*(z/z0)^0.40
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassF},::Type{CCPSUrban}) = u0*(z/z0)^0.60


"""
    crosswind_dispersion(x, Plume, StabilityClass, CCPSUrban)

Plume crosswind dispersion correlations, for urban terrain

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
+ Briggs, Gary A. 1973. *Diffusion Estimation for Small Emissions. Preliminary Report*. United States. https://doi.org/10.2172/5118833
+ Griffiths, R. F. 1994. "Errors in the use of the Briggs parameterization for atmospheric dispersion coefficients." *Atmospheric Environment* 28(17):2861-2865 https://doi.org/10.1016/1352-2310(94)90086-8
"""
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::Type{CCPSUrban}) = 0.32x/√(1+0.0004x)
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::Type{CCPSUrban}) = 0.32x/√(1+0.0004x)
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::Type{CCPSUrban}) = 0.22x/√(1+0.0004x)
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::Type{CCPSUrban}) = 0.16x/√(1+0.0004x)
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::Type{CCPSUrban}) = 0.11x/√(1+0.0004x)
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::Type{CCPSUrban}) = 0.11x/√(1+0.0004x)


"""
    vertical_dispersion(x, Plume, StabilityClass, CCPSUrban)

Plume vertical dispersion correlations, for urban terrain

References:
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
+ Briggs, Gary A. 1973. *Diffusion Estimation for Small Emissions. Preliminary Report*. United States. https://doi.org/10.2172/5118833
+ Griffiths, R. F. 1994. "Errors in the use of the Briggs parameterization for atmospheric dispersion coefficients." *Atmospheric Environment* 28(17):2861-2865 https://doi.org/10.1016/1352-2310(94)90086-8
"""
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::Type{CCPSUrban}) = 0.24x*√(1+0.001x)
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::Type{CCPSUrban}) = 0.24x*√(1+0.001x)
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::Type{CCPSUrban}) = 0.20x
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::Type{CCPSUrban}) = 0.14x/√(1+0.0003x)
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::Type{CCPSUrban}) = 0.08x/√(1+0.0015x)
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::Type{CCPSUrban}) = 0.08x/√(1+0.0015x)


# Puff dispersion correlations
"""
    crosswind_dispersion(x, Puff, StabilityClass. {CCPSUrban,CCPSRural})

Puff crosswind dispersion correlations

References:
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
"""
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassA}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.18*x^0.92
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassB}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.14*x^0.92
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassC}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.10*x^0.92
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassD}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.06*x^0.92
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassE}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.04*x^0.92
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassF}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.02*x^0.89


"""
    vertical_dispersion(x, Puff, StabilityClass, {CCPSUrban,CCPSRural})

Puff vertical dispersion correlations

References:
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
"""
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassA}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.60*x^0.75
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassB}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.53*x^0.73
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassC}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.34*x^0.71
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassD}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.15*x^0.70
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassE}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.10*x^0.65
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassF}, ::Union{Type{CCPSUrban},Type{CCPSRural}}) = 0.05*x^0.61
