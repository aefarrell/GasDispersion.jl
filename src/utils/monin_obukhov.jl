"""
    _monin_obukhov(roughness, StabilityClass, WindEquation)
returns the Monin-Obukhov length for a given Pasquill-Gifford stability class
and surface roughness (in meters)

# References
+ Pasquill, Frank. 1974. *Atmospheric Diffusion, 2nd Ed.* New York: Halstead Press, New York

"""
_monin_obukhov(zR::Number, ::Type{ClassA}, ::Any) = -11.4*zR^0.10
_monin_obukhov(zR::Number, ::Type{ClassB}, ::Any) = -26.0*zR^0.17
_monin_obukhov(zR::Number, ::Type{ClassC}, ::Any) = -123.0*zR^0.30
_monin_obukhov(zR::Number, ::Type{ClassD}, ::Any) = Inf
_monin_obukhov(zR::Number, ::Type{ClassE}, ::Any) = 123.0*zR^0.30
_monin_obukhov(zR::Number, ::Type{ClassF}, ::Any) = 26.0*zR^0.17

_monin_obukhov(zR::Number, s::Any) = _monin_obukhov(zR, s, Nothing)