# Correlations from the CCPS *Guidelines for Consequence Analysis of Chemical Releases*


# Correlations for rural terrain
CCPSRural() = BasicEquationSet(IrwinRural(),nothing,BriggsRuralσy(),BriggsRuralσz())

# Correlations for urban terrain
CCPSUrban() = BasicEquationSet(IrwinUrban(),nothing,BriggsUrbanσy(),BriggsUrbanσz())

# Correlations for puffs, rural terrain
CCPSPuffRural() = BasicEquationSet(IrwinRural(),CCPSPuffσx(),CCPSPuffσy(),CCPSPuffσz())

# Correlations for puffs, urban terrain
CCPSPuffUrban() = BasicEquationSet(IrwinUrban(),CCPSPuffσx(),CCPSPuffσy(),CCPSPuffσz())