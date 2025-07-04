# Correlations for rural terrain
struct IrwinRural <: PowerLawWind end

"""
    windspeed(u0::Number,z0::Number,z::Number, stability, CCPSRural)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

Assumes rural terrain.

# References
+ Irwin, J. S. 1979. "A Theoretical Variation of the Wind Profile Power Law Exponent as a Function of Surface Roughness, and Stability," *Atmospheric Environment, 13: 191-194.
"""
windspeed(u0::Number,z0::Number,z::Number,::ClassA,::IrwinRural) = u0*(z/z0)^0.07
windspeed(u0::Number,z0::Number,z::Number,::ClassB,::IrwinRural) = u0*(z/z0)^0.07
windspeed(u0::Number,z0::Number,z::Number,::ClassC,::IrwinRural) = u0*(z/z0)^0.10
windspeed(u0::Number,z0::Number,z::Number,::ClassD,::IrwinRural) = u0*(z/z0)^0.15
windspeed(u0::Number,z0::Number,z::Number,::ClassE,::IrwinRural) = u0*(z/z0)^0.35
windspeed(u0::Number,z0::Number,z::Number,::ClassF,::IrwinRural) = u0*(z/z0)^0.55


# Correlations for urban terrain
struct IrwinUrban <: PowerLawWind end

"""
    windspeed(u0::Number,z0::Number,z::Number, stability, CCPSUrban)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

Assumes urban terrain.

# References
+ Irwin, J. S. 1979. "A Theoretical Variation of the Wind Profile Power Law Exponent as a Function of Surface Roughness, and Stability," *Atmospheric Environment, 13: 191-194.
"""
windspeed(u0::Number,z0::Number,z::Number,::ClassA,::IrwinUrban) = u0*(z/z0)^0.15
windspeed(u0::Number,z0::Number,z::Number,::ClassB,::IrwinUrban) = u0*(z/z0)^0.15
windspeed(u0::Number,z0::Number,z::Number,::ClassC,::IrwinUrban) = u0*(z/z0)^0.20
windspeed(u0::Number,z0::Number,z::Number,::ClassD,::IrwinUrban) = u0*(z/z0)^0.25
windspeed(u0::Number,z0::Number,z::Number,::ClassE,::IrwinUrban) = u0*(z/z0)^0.40
windspeed(u0::Number,z0::Number,z::Number,::ClassF,::IrwinUrban) = u0*(z/z0)^0.60
