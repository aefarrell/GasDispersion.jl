abstract type PropertyCorrelation end
abstract type EquationNumber end

struct Eq1 <: EquationNumber end
struct Eq2 <: EquationNumber end

struct Antoine{F<:Number} <: PropertyCorrelation
    A::F
    B::F
    C::F
end
Antoine(A,B,C) = Antoine(promote(A,B,C)...)

function (a::Antoine)(T)
    return exp(a.A - a.B/(T + a.C))
end

struct DIPPRVaporPressure{F<:Number} <: PropertyCorrelation
    C1::F
    C2::F
    C3::F
    C4::F
    C5::F
end
DIPPRVaporPressure(C1,C2,C3,C4,C5) = DIPPRVaporPressure(promote(C1,C2,C3,C4,C5)...)

function (p::DIPPRVaporPressure)(T)
    lnP = p.C1 + p.C2/T + p.C3*log(T) + p.C4*T^p.C5
    return exp(lnP)
end

struct DIPPRLiquidDensity{F<:Number} <: PropertyCorrelation
    MW::F
    Tc::F
    C1::F
    C2::F
    C3::F
    C4::F
end
DIPPRLiquidDensity(MW,Tc,C1,C2,C3,C4) = DIPPRLiquidDensity(promote(MW,Tc,C1,C2,C3,C4)...)

function (l::DIPPRLiquidDensity)(T,P)
    #Tr = T/l.Tc
    ρl = l.C1/( l.C2^(1 + (1 - T/l.C3)^l.C4) )
    return ρl*l.MW*1000
end

struct DIPPRLatentHeat{F<:Number} <: PropertyCorrelation
    MW::F
    Tc::F
    C1::F
    C2::F
    C3::F
    C4::F
    C5::F
end
DIPPRLatentHeat(MW,Tc,C1,C2,C3,C4,C5) = DIPPRLatentHeat(promote(MW,Tc,C1,C2,C3,C4,C5)...)

function (l::DIPPRLatentHeat)(T)
    Tr = T/l.Tc
    ΔHv = l.C1*(1-Tr)^(l.C2 + l.C3*Tr + l.C4*Tr^2 + l.C5*Tr^3)
    return ΔHv/(l.MW*1000)
end

struct DIPPRIdealGasHeatCapacity{F<:Number} <: PropertyCorrelation
    MW::F
    Tc::F
    C1::F
    C2::F
    C3::F
    C4::F
    C5::F
end
DIPPRIdealGasHeatCapacity(MW,Tc,C1,C2,C3,C4,C5) = DIPPRIdealGasHeatCapacity(promote(MW,Tc,C1,C2,C3,C4,C5)...)

function (c::DIPPRIdealGasHeatCapacity)(T)
    #Tr = T/c.Tc
    Cp = c.C1 + c.C2*((c.C3/T)/sinh(c.C3/T))^2 + c.C4*((c.C5/T)/cosh(c.C5/T))^2
    return Cp/(c.MW*1000)
end

struct DIPPRLiquidHeatCapacity{E<:EquationNumber, F<:Number} <: PropertyCorrelation
    EQ::Type{E}
    MW::F
    Tc::F
    C1::F
    C2::F
    C3::F
    C4::F
    C5::F
end
DIPPRLiquidHeatCapacity(EQ,MW,Tc,C1,C2,C3,C4,C5) = DIPPRLiquidHeatCapacity(EQ,promote(MW,Tc,C1,C2,C3,C4,C5)...)

function (c::DIPPRLiquidHeatCapacity{Eq1,<:Number})(T)
    Cp = c.C1 + c.C2*T + c.C3*T^2 + c.C4*T^3 + c.C5*T^4
    return Cp/(c.MW*1000)
end

function (c::DIPPRLiquidHeatCapacity{Eq2,<:Number})(T)
    Tr = T/c.Tc
    t = 1-Tr
    Cp = c.C1^2/t + c.C2 - 2*c.C1*c.C3*t - c.C1*c.C4*t^2 - c.C3^2*t^3/3 - c.C3*c.C4*t^4/2 - c.C4^2*t^5/5
    return Cp/(c.MW*1000)
end