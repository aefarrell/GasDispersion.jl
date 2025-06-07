"""
    monin_obuknov(roughness, StabilityClass, WindEquation)
returns the Monin-Obukhov length for a given Pasquill-Gifford stability class
and surface roughness (in meters)

# References
+ Pasquill, Frank. 1974. *Atmospheric Diffusion, 2nd Ed.* New York: Halstead Press, New York

"""
monin_obuknov(zR::Number, ::Type{ClassA}, ::Any) = -11.4*zR^0.10
monin_obuknov(zR::Number, ::Type{ClassB}, ::Any) = -26.0*zR^0.17
monin_obuknov(zR::Number, ::Type{ClassC}, ::Any) = -123.0*zR^0.30
monin_obuknov(zR::Number, ::Type{ClassD}, ::Any) = Inf
monin_obuknov(zR::Number, ::Type{ClassE}, ::Any) = 123.0*zR^0.30
monin_obuknov(zR::Number, ::Type{ClassF}, ::Any) = 26.0*zR^0.17

monin_obuknov(zR::Number, s::Any) = monin_obuknov(zR, s, Nothing)