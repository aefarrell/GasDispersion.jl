"""
    _monin_obukhov(roughness, StabilityClass)
returns the Monin-Obukhov length for a given Pasquill-Gifford stability class
and surface roughness (in meters)

# References
+ Pasquill, F., *Atmospheric Diffusion, 2nd Ed.*, Halstead Press, New York (1974)

"""
_monin_obukhov(zR::Number, ::Type{ClassA}) = -11.4*zR^0.10
_monin_obukhov(zR::Number, ::Type{ClassB}) = -26.0*zR^0.17
_monin_obukhov(zR::Number, ::Type{ClassC}) = -123.0*zR^0.30
_monin_obukhov(zR::Number, ::Type{ClassD}) = Inf
_monin_obukhov(zR::Number, ::Type{ClassE}) = 123.0*zR^0.30
_monin_obukhov(zR::Number, ::Type{ClassF}) = 26.0*zR^0.17
