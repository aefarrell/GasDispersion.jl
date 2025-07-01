# correlations for Monin-Obukhov length
include("monin_obukhov.jl")

# correlations for windspeed
include("wind_profile.jl")

# Pasquill-Gifford dispersion correlations
include("pasquill_gifford.jl")

# Briggs' model for plume rise
include("plume_rise.jl")

# mixing height correlations
include("mixing_layer.jl")

# model correlations
include("britter_mcquaid_correls.jl")

## Equation Sets
# includes for defined sets of correlations go here
include("equation_sets/briggs.jl")
include("equation_sets/irwin.jl")
include("equation_sets/ccps.jl")
include("equation_sets/tno.jl")
include("equation_sets/turner.jl")
include("equation_sets/isc3.jl")

# Default correlations
DefaultSet = BasicEquationSet{DefaultWind,Nothing,Defaultσy,Defaultσz}
DefaultPuffSet = BasicEquationSet{DefaultWind,CCPSPuffσx,CCPSPuffσy,CCPSPuffσz}

# a lazy check for if x,y,z are in the domain
struct ProblemDomain{F<:Number}
    xmin::F
    xmax::F
    ymin::F
    ymax::F
    zmin::F
    zmax::F
end

function _in_domain(x::Number, y::Number, z::Number, domain::ProblemDomain)
    return x ≥ domain.xmin && x ≤ domain.xmax &&
           y ≥ domain.ymin && y ≤ domain.ymax &&
           z ≥ domain.zmin && z ≤ domain.zmax
end