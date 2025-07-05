# Correlations from the CCPS *Guidelines for Consequence Analysis of Chemical Releases*


# Correlations for rural terrain
const CCPSRural = BasicEquationSet(IrwinRural(),nothing,BriggsRuralσy(),BriggsRuralσz())

# Correlations for urban terrain
const CCPSUrban = BasicEquationSet(IrwinUrban(),nothing,BriggsUrbanσy(),BriggsUrbanσz())

# Correlations for puffs, rural terrain
const CCPSPuffRural = BasicEquationSet(IrwinRural(),CCPSPuffσx(),CCPSPuffσy(),CCPSPuffσz())

# Correlations for puffs, urban terrain
const CCPSPuffUrban = BasicEquationSet(IrwinUrban(),CCPSPuffσx(),CCPSPuffσy(),CCPSPuffσz())