# Correlations from the CCPS *Guidelines for Consequence Analysis of Chemical Releases*


# Correlations for rural terrain
CCPSRural = BasicEquationSet{IrwinRural,Nothing,BriggsRuralσy,BriggsRuralσz}

# Correlations for urban terrain
CCPSUrban = BasicEquationSet{IrwinUrban,Nothing,BriggsUrbanσy,BriggsUrbanσz}

# Correlations for puffs, rural terrain
CCPSPuffRural = BasicEquationSet{IrwinRural,CCPSPuffσx,CCPSPuffσy,CCPSPuffσz}

# Correlations for puffs, urban terrain
CCPSPuffUrban = BasicEquationSet{IrwinUrban,CCPSPuffσx,CCPSPuffσy,CCPSPuffσz}