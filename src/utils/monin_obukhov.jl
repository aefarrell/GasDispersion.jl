"""
    _monin_obukhov(roughness, StabilityClass)
returns the Monin-Obukhov length for a given Pasquill-Gifford stability class
and surface roughness (in meters)
Curve fit from
    Pasquill, F., *Atmospheric Diffusion, 2nd Ed.*, Halstead Press, New York, 1974.
"""
_monin_obukhov(zR, ::Type{ClassA}) = -11.4*zR^0.10
_monin_obukhov(zR, ::Type{ClassB}) = -26.0*zR^0.17
_monin_obukhov(zR, ::Type{ClassC}) = -123.0*zR^0.30
_monin_obukhov(zR, ::Type{ClassD}) = Inf
_monin_obukhov(zR, ::Type{ClassE}) = 123.0*zR^0.30
_monin_obukhov(zR, ::Type{ClassF}) = 26.0*zR^0.17
