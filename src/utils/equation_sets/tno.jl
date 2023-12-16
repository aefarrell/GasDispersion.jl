# Correlations from the TNO Yellow Book

struct TNO <: EquationSet end

"""
    crosswind_dispersion(x, Plume, StabilityClass, TNO)

Plume crosswind dispersion correlations

# References
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::Type{TNO}) = 0.527x^0.865
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::Type{TNO}) = 0.371x^0.866
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::Type{TNO}) = 0.209x^0.897
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::Type{TNO}) = 0.128x^0.905
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::Type{TNO}) = 0.098x^0.902
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::Type{TNO}) = 0.065x^0.902


"""
    vertical_dispersion(x, Plume, StabilityClass, TNO)

Plume vertical dispersion correlations

References:
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::Type{TNO}) = 0.28x^0.90
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::Type{TNO}) = 0.23x^0.85
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::Type{TNO}) = 0.22x^0.80
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::Type{TNO}) = 0.20x^0.76
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::Type{TNO}) = 0.15x^0.73
vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::Type{TNO}) = 0.12x^0.67

"""
    crosswind_dispersion(x, Puff, StabilityClass, TNO)

Plume crosswind dispersion correlations

# References
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassA}, es::Type{TNO}) = 0.5*crosswind_dispersion(x,Plume,ClassA,es)
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassB}, es::Type{TNO}) = 0.5*crosswind_dispersion(x,Plume,ClassB,es)
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassC}, es::Type{TNO}) = 0.5*crosswind_dispersion(x,Plume,ClassC,es)
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassD}, es::Type{TNO}) = 0.5*crosswind_dispersion(x,Plume,ClassD,es)
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassE}, es::Type{TNO}) = 0.5*crosswind_dispersion(x,Plume,ClassE,es)
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassF}, es::Type{TNO}) = 0.5*crosswind_dispersion(x,Plume,ClassF,es)

"""
    vertical_dispersion(x, Puff, StabilityClass, TNO)

Puff vertical dispersion correlations

References:
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassA}, ::Type{TNO}) = 0.28x
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassB}, ::Type{TNO}) = 0.23x
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassC}, ::Type{TNO}) = 0.22x
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassD}, ::Type{TNO}) = 0.20x
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassE}, ::Type{TNO}) = 0.15x
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassF}, ::Type{TNO}) = 0.12x

"""
    downwind_dispersion(x, Puff, StabilityClass, TNO)

Puff downwind dispersion correlations

References:
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
"""
downwind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassA}, ::Type{TNO}) = 0.13x
downwind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassB}, ::Type{TNO}) = 0.13x
downwind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassC}, ::Type{TNO}) = 0.13x
downwind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassD}, ::Type{TNO}) = 0.13x
downwind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassE}, ::Type{TNO}) = 0.13x
downwind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassF}, ::Type{TNO}) = 0.13x