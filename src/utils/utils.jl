# correlations for Monin-Obukhov length
include("monin_obukhov.jl")

# correlations for windspeed
include("wind_profile.jl")

# Pasquill-Gifford dispersion correlations
include("pasquill_gifford.jl")

# Briggs' model for plume rise
include("plume_rise.jl")

# model correlations
include("britter_mcquaid_correls.jl")

## Equation Sets
# includes for defined sets of correlations go here
include("equation_sets/briggs.jl")
include("equation_sets/irwin.jl")
include("equation_sets/ccps.jl")
# include("equation_sets/tno.jl")
include("equation_sets/turner.jl")
include("equation_sets/isc3.jl")

# Default correlations
DefaultSet = BasicEquationSet{DefaultWind,Nothing,Defaultσy,Defaultσz}
DefaultPuffSet = BasicEquationSet{DefaultWind,CCPSPuffσx,CCPSPuffσy,CCPSPuffσz}