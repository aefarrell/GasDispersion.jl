# Correlations from the TNO Yellow Book

struct TNO <: EquationSet end

# no power law correlation used, just passes to default
_windspeed(a::Atmosphere,z::Number,::TNO) = _windspeed(a,z,DefaultSet())
_windspeed(u0::Number,z0::Number,z::Number,s::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD},Type{ClassE},Type{ClassF}}, ::TNO) = _windspeed(u0,z0,z,s,DefaultSet())


"""
    crosswind_dispersion(x, Plume, StabilityClass, ::TNO)

Plume crosswind dispersion correlations

# References
+ TNO Yellow Book, section 4.5.3.4 (1)
"""
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::TNO) = 0.527x^0.865
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::TNO) = 0.371x^0.866
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::TNO) = 0.209x^0.897
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::TNO) = 0.128x^0.905
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::TNO) = 0.098x^0.902
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::TNO) = 0.065x^0.902

"""
    crosswind_dispersion(x, Puff, StabilityClass, ::TNO)

Puff crosswind dispersion correlations

# References
+ TNO Yellow Book, section 4.5.3.4 (4)
"""
crosswind_dispersion(x::Number, ::Type{Puff}, s::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD},Type{ClassE},Type{ClassF}}, eq::TNO) = 0.5*crosswind_dispersion(x, Plume, s, eq)


"""
    vertical_dispersion(x, Plume, StabilityClass, ::TNO)

Plume vertical dispersion correlations

References:
+ TNO Yellow Book, section 4.5.3.4 (1)
"""
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::TNO) = 0.28x^0.90
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::TNO) = 0.23x^0.85
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::TNO) = 0.22x^0.80
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::TNO) = 0.20x^0.76
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::TNO) = 0.15x^0.73
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::TNO) = 0.12x^0.67

"""
    vertical_dispersion(x, Puff, StabilityClass, ::TNO)

Puff vertical dispersion correlations

References:
+ TNO Yellow Book, section 4.5.3.4 (1)
"""
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassA}, ::TNO) = 0.28x
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassB}, ::TNO) = 0.23x
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassC}, ::TNO) = 0.22x
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassD}, ::TNO) = 0.20x
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassE}, ::TNO) = 0.15x
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassF}, ::TNO) = 0.12x

"""
    downwind_dispersion(x, Puff, StabilityClass)

Puff downwind dispersion correlations

References:
+ TNO Yellow Book, section 4.5.3.4 (1)
"""
downwind_dispersion(x::Number, ::Type{Puff}, ::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD},Type{ClassE},Type{ClassF}}, ::TNO) = 0.13x