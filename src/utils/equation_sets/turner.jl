# Correlations from Turner (1970)

struct Turner <: EquationSet end

"""
    crosswind_dispersion(x, Plume, StabilityClass, Turner)

Plume crosswind dispersion correlations

# References
+ Lees, Frank P. 1996. *Loss Prevention in the Process Industries, 2nd ed*. Oxford: Butterworth-Heinemann
+ Turner, D. Bruce. 1970. *Workbook of Atmospheric Dispersion Estimates*. United States Environmental Protection Agency.
"""
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::Type{Turner}) = 0.493x^0.88
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::Type{Turner}) = 0.337x^0.88
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::Type{Turner}) = 0.195x^0.90
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::Type{Turner}) = 0.128x^0.90
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::Type{Turner}) = 0.091x^0.91
crosswind_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::Type{Turner}) = 0.067x^0.90

"""
    vertical_dispersion(x, Plume, StabilityClass, Turner)

Plume vertical dispersion correlations

References:
+ Lees, Frank P. 1996. *Loss Prevention in the Process Industries, 2nd ed*. Oxford: Butterworth-Heinemann
+ Turner, D. Bruce. 1970. *Workbook of Atmospheric Dispersion Estimates*. United States Environmental Protection Agency.
"""
function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassA}, ::Type{Turner})
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

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassB}, ::Type{Turner})
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

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassC}, ::Type{Turner})
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

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassD}, ::Type{Turner})
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

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassE}, ::Type{Turner})
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

function vertical_dispersion(x::Number, ::Type{Plume}, ::Type{ClassF}, ::Type{Turner})
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
