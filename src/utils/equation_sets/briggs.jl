# Correlations from Briggs, Gary A. 1973. Diffusion Estimation for Small Emissions. Preliminary Report. United States


# Correlations for rural terrain
struct BriggsRuralσy <: DispersionFunction end
struct BriggsRuralσz <: DispersionFunction end

"""
    crosswind_dispersion(x, StabilityClass, Briggsσy)

Plume crosswind dispersion correlations, for rural terrain

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
+ Briggs, Gary A. 1973. *Diffusion Estimation for Small Emissions. Preliminary Report*. United States. https://doi.org/10.2172/5118833
+ Griffiths, R. F. 1994. "Errors in the use of the Briggs parameterization for atmospheric dispersion coefficients." *Atmospheric Environment* 28(17):2861-2865 https://doi.org/10.1016/1352-2310(94)90086-8
"""
crosswind_dispersion(x::Number, ::ClassA, ::BriggsRuralσy) = 0.22x/√(1+0.0001x)
crosswind_dispersion(x::Number, ::ClassB, ::BriggsRuralσy) = 0.16x/√(1+0.0001x)
crosswind_dispersion(x::Number, ::ClassC, ::BriggsRuralσy) = 0.11x/√(1+0.0001x)
crosswind_dispersion(x::Number, ::ClassD, ::BriggsRuralσy) = 0.08x/√(1+0.0001x)
crosswind_dispersion(x::Number, ::ClassE, ::BriggsRuralσy) = 0.06x/√(1+0.0001x)
crosswind_dispersion(x::Number, ::ClassF, ::BriggsRuralσy) = 0.04x/√(1+0.0001x)


"""
    vertical_dispersion(x, StabilityClass, BriggsRuralσz)

Plume vertical dispersion correlations, for rural terrain

References:
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
+ Briggs, Gary A. 1973. *Diffusion Estimation for Small Emissions. Preliminary Report*. United States. https://doi.org/10.2172/5118833
+ Griffiths, R. F. 1994. "Errors in the use of the Briggs parameterization for atmospheric dispersion coefficients." *Atmospheric Environment* 28(17):2861-2865 https://doi.org/10.1016/1352-2310(94)90086-8
"""
vertical_dispersion(x::Number, ::ClassA, ::BriggsRuralσz) = 0.20x
vertical_dispersion(x::Number, ::ClassB, ::BriggsRuralσz) = 0.12x
vertical_dispersion(x::Number, ::ClassC, ::BriggsRuralσz) = 0.08x/√(1+0.0002x)
vertical_dispersion(x::Number, ::ClassD, ::BriggsRuralσz) = 0.06x/√(1+0.0015x)
vertical_dispersion(x::Number, ::ClassE, ::BriggsRuralσz) = 0.03x/(1+0.0003x)
vertical_dispersion(x::Number, ::ClassF, ::BriggsRuralσz) = 0.016x/(1+0.0003x)


# Correlations for urban terrain
struct BriggsUrbanσy <: DispersionFunction end
struct BriggsUrbanσz <: DispersionFunction end

"""
    crosswind_dispersion(x, StabilityClass, BriggsUrbanσy)

Plume crosswind dispersion correlations, for urban terrain

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
+ Briggs, Gary A. 1973. *Diffusion Estimation for Small Emissions. Preliminary Report*. United States. https://doi.org/10.2172/5118833
+ Griffiths, R. F. 1994. "Errors in the use of the Briggs parameterization for atmospheric dispersion coefficients." *Atmospheric Environment* 28(17):2861-2865 https://doi.org/10.1016/1352-2310(94)90086-8
"""
crosswind_dispersion(x::Number, ::ClassA, ::BriggsUrbanσy) = 0.32x/√(1+0.0004x)
crosswind_dispersion(x::Number, ::ClassB, ::BriggsUrbanσy) = 0.32x/√(1+0.0004x)
crosswind_dispersion(x::Number, ::ClassC, ::BriggsUrbanσy) = 0.22x/√(1+0.0004x)
crosswind_dispersion(x::Number, ::ClassD, ::BriggsUrbanσy) = 0.16x/√(1+0.0004x)
crosswind_dispersion(x::Number, ::ClassE, ::BriggsUrbanσy) = 0.11x/√(1+0.0004x)
crosswind_dispersion(x::Number, ::ClassF, ::BriggsUrbanσy) = 0.11x/√(1+0.0004x)


"""
    vertical_dispersion(x, StabilityClass, BriggsUrbanσz)

Plume vertical dispersion correlations, for urban terrain

References:
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
+ Briggs, Gary A. 1973. *Diffusion Estimation for Small Emissions. Preliminary Report*. United States. https://doi.org/10.2172/5118833
+ Griffiths, R. F. 1994. "Errors in the use of the Briggs parameterization for atmospheric dispersion coefficients." *Atmospheric Environment* 28(17):2861-2865 https://doi.org/10.1016/1352-2310(94)90086-8
"""
vertical_dispersion(x::Number, ::ClassA, ::BriggsUrbanσz) = 0.24x*√(1+0.001x)
vertical_dispersion(x::Number, ::ClassB, ::BriggsUrbanσz) = 0.24x*√(1+0.001x)
vertical_dispersion(x::Number, ::ClassC, ::BriggsUrbanσz) = 0.20x
vertical_dispersion(x::Number, ::ClassD, ::BriggsUrbanσz) = 0.14x/√(1+0.0003x)
vertical_dispersion(x::Number, ::ClassE, ::BriggsUrbanσz) = 0.08x/√(1+0.0015x)
vertical_dispersion(x::Number, ::ClassF, ::BriggsUrbanσz) = 0.08x/√(1+0.0015x)