# Correlations for rural terrain
struct IrwinRural <: PowerLawWind end

"""
    _windspeed(u0::Number,z0::Number,z::Number, stability, CCPSRural)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

Assumes rural terrain.

# References
+ Irwin, J. S. 1979. "A Theoretical Variation of the Wind Profile Power Law Exponent as a Function of Surface Roughness, and Stability," *Atmospheric Environment, 13: 191-194.
"""
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassA},::Type{IrwinRural}) = u0*(z/z0)^0.07
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassB},::Type{IrwinRural}) = u0*(z/z0)^0.07
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassC},::Type{IrwinRural}) = u0*(z/z0)^0.10
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassD},::Type{IrwinRural}) = u0*(z/z0)^0.15
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassE},::Type{IrwinRural}) = u0*(z/z0)^0.35
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassF},::Type{IrwinRural}) = u0*(z/z0)^0.55


# Correlations for urban terrain
struct IrwinUrban <: PowerLawWind end

"""
    _windspeed(u0::Number,z0::Number,z::Number, stability, CCPSUrban)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

Assumes urban terrain.

# References
+ Irwin, J. S. 1979. "A Theoretical Variation of the Wind Profile Power Law Exponent as a Function of Surface Roughness, and Stability," *Atmospheric Environment, 13: 191-194.
"""
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassA},::Type{IrwinUrban}) = u0*(z/z0)^0.15
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassB},::Type{IrwinUrban}) = u0*(z/z0)^0.15
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassC},::Type{IrwinUrban}) = u0*(z/z0)^0.20
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassD},::Type{IrwinUrban}) = u0*(z/z0)^0.25
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassE},::Type{IrwinUrban}) = u0*(z/z0)^0.40
_windspeed(u0::Number,z0::Number,z::Number,::Type{ClassF},::Type{IrwinUrban}) = u0*(z/z0)^0.60
