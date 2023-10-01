"""
    crosswind_dispersion(x, Plume, StabilityClass; avg_time=600.0)

Plume crosswind dispersion correlations

# References
+ Spicer, Thomas O. and Jerry A. Havens. 1988. *Development of Vapor Dispersion Models for Non-Neutrally Buoyant Gas Mixtures--Analysis of TFI/NH3 Test Data*. United States.

"""
function crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::DefaultSet; avg_time=600.0)
    δ, β, tₐ = 0.423, 0.9, 18.4
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end

function crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::DefaultSet; avg_time=600.0)
    δ, β, tₐ = 0.313, 0.9, 18.4
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end

function crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::DefaultSet; avg_time=600.0)
    δ, β, tₐ = 0.210, 0.9, 18.4
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end

function crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::DefaultSet; avg_time=600.0)
    δ, β, tₐ = 0.136, 0.9, 18.3
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end

function crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::DefaultSet; avg_time=600.0)
    δ, β, tₐ = 0.102, 0.9, 11.4
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end

function crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::DefaultSet; avg_time=600.0)
    δ, β, tₐ = 0.0674, 0.9, 4.6
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end


"""
    vertical_dispersion(x, Plume, StabilityClass)

Plume vertical dispersion correlations

References:
+ Seinfeld, John H. 1986. *Atmospheric Chemistry and Physics of Air Pollution*. New York: John Wiley and Sons
"""
function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::DefaultSet)
    δ = 107.7
    β = -1.7172
    γ = 0.2770
    return δ*(x^β)*exp(γ*log(x)^2)
end

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::DefaultSet)
    δ=0.1355
    β=0.8752
    γ=0.0136
    return δ*(x^β)*exp(γ*log(x)^2)
end

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::DefaultSet)
    δ=0.09623
    β=0.9477
    γ=-0.0020
    return δ*(x^β)*exp(γ*log(x)^2)
end

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::DefaultSet)
    δ=0.04134
    β=1.1737
    γ=-0.0316
    return δ*(x^β)*exp(γ*log(x)^2)
end

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::DefaultSet)
    δ=0.02275
    β=1.3010
    γ=-0.0450
    return δ*(x^β)*exp(γ*log(x)^2)
end

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::DefaultSet)
    δ=0.01122
    β=1.4024
    γ=-0.0540
    return δ*(x^β)*exp(γ*log(x)^2)
end

"""
    crosswind_dispersion(x, Puff, StabilityClass)

Puff crosswind dispersion correlations

References:
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
"""
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassA}, ::DefaultSet) = 0.18*x^0.92
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassB}, ::DefaultSet) = 0.14*x^0.92
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassC}, ::DefaultSet) = 0.10*x^0.92
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassD}, ::DefaultSet) = 0.06*x^0.92
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassE}, ::DefaultSet) = 0.04*x^0.92
crosswind_dispersion(x::Number, ::Type{Puff}, ::Type{ClassF}, ::DefaultSet) = 0.02*x^0.89


"""
    downwind_dispersion(x, Puff, StabilityClass)

Puff downwind dispersion correlations

References:
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
"""
function downwind_dispersion(x::Number, ::Type{Puff}, stab::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD},Type{ClassE},Type{ClassF}}, ::DefaultSet)
    return crosswind_dispersion(x, Puff, stab, DefaultSet())
end


"""
    vertical_dispersion(x, Puff, StabilityClass)

Puff vertical dispersion correlations

References:
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
"""
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassA}, ::DefaultSet) = 0.60*x^0.75
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassB}, ::DefaultSet) = 0.53*x^0.73
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassC}, ::DefaultSet) = 0.34*x^0.71
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassD}, ::DefaultSet) = 0.15*x^0.70
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassE}, ::DefaultSet) = 0.10*x^0.65
vertical_dispersion(x::Number, ::Type{Puff}, ::Type{ClassF}, ::DefaultSet) = 0.05*x^0.61
