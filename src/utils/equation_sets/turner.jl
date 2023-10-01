# Correlations from Turner (1970)

struct Turner <: EquationSet end

# no power law correlation used, just passes to default
_windspeed(a::Atmosphere,z::Number,::Turner) = _windspeed(a,z,DefaultSet())
_windspeed(u0::Number,z0::Number,z::Number,s::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD},Type{ClassE},Type{ClassF}},::Turner) = _windspeed(u0,z0,z,s,DefaultSet())

# no puff correlations used, just passes to default
crosswind_dispersion(x::Number,::Type{Puff},s::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD},Type{ClassE},Type{ClassF}},::Turner) = crosswind_dispersion(x,Puff,s,DefaultSet())
vertical_dispersion(x::Number,::Type{Puff},s::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD},Type{ClassE},Type{ClassF}},::Turner) = vertical_dispersion(x,Puff,s,DefaultSet())
downwind_dispersion(x::Number,::Type{Puff},s::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD},Type{ClassE},Type{ClassF}},::Turner) = downwind_dispersion(x,Puff,s,DefaultSet())


"""
    crosswind_dispersion(x, Plume, StabilityClass, ::Turner)

Plume crosswind dispersion correlations

# References
+ Lees, Frank P. 1996. *Loss Prevention in the Process Industries, 2nd ed*. Oxford: Butterworth-Heinemann
+ Turner, D. Bruce. 1970. *Workbook of Atmospheric Dispersion Estimates*. United States Environmental Protection Agency.
"""
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::Turner) = 0.493x^0.88
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::Turner) = 0.337x^0.88
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::Turner) = 0.195x^0.90
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::Turner) = 0.128x^0.90
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::Turner) = 0.091x^0.91
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::Turner) = 0.067x^0.90

"""
    vertical_dispersion(x, Plume, StabilityClass, ::Turner)

Plume vertical dispersion correlations

References:
+ Lees, Frank P. 1996. *Loss Prevention in the Process Industries, 2nd ed*. Oxford: Butterworth-Heinemann
+ Turner, D. Bruce. 1970. *Workbook of Atmospheric Dispersion Estimates*. United States Environmental Protection Agency.
"""
function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::Turner)
    if x < 100
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 3000"
        return 0.087x^1.10
    elseif x ≤ 300
        return 0.087x^1.10
    elseif x ≤ 3000
        logsz = -1.67+0.902log10(x)+0.181(log10(x))^2
        return 10^logsz
    else
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 3000"
        logsz = -1.67+0.902log10(x)+0.181(log10(x))^2
        return 10^logsz
    end
end

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::Turner)
    if x < 100
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 20000"
        return 0.135x^0.95
    elseif x ≤ 500
        return 0.135x^0.95
    elseif x ≤ 20000
        logsz = -1.25+1.09log10(x)+0.0018(log10(x))^2
        return 10^logsz
    else
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 20000"
        logsz = -1.25+1.09log10(x)+0.0018(log10(x))^2
        return 10^logsz
    end
end

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::Turner)
    if x < 100
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 100000"
        return 0.112x^0.91
    elseif x ≤ 100000
        return 0.112x^0.91
    else
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 100000"
        return 0.112x^0.91
    end
end

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::Turner)
    if x < 100
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 100000"
        return 0.093x^0.85
    elseif x ≤ 500
        return 0.093x^0.85
    elseif x ≤ 100000
        logsz = -1.22+1.08log10(x)-0.061(log10(x))^2
        return 10^logsz
    else
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 100000"
        logsz = -1.22+1.08log10(x)-0.061(log10(x))^2
        return 10^logsz
    end
end

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::Turner)
    if x < 100
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 100000"
        return 0.082x^0.82
    elseif x ≤ 500
        return 0.082x^0.82
    elseif x ≤ 100000
        logsz = -1.18+1.04log10(x)-0.070(log10(x))^2
        return 10^logsz
    else
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 100000"
        logsz = -1.18+1.04log10(x)-0.070(log10(x))^2
        return 10^logsz
    end
end

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::Turner)
    if x < 100
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 100000"
        return 0.057x^0.80
    elseif x ≤ 500
        return 0.057x^0.80
    elseif x ≤ 100000
        logsz = -1.91+1.37log10(x)-0.119(log10(x))^2
        return 10^logsz
    else
        @warn "x = $x, outside of the range of the Turner correlation, range is 100 ≤ x ≤ 100000"
        logsz = -1.91+1.37log10(x)-0.119(log10(x))^2
        return 10^logsz
    end
end
