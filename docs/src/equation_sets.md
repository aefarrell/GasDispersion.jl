# Equation Sets

The various dispersion models each depend upon several parameters which are themselves, 
often, correlations. For any given parameter there are several different correlations 
in the literature. To make this more transparent, sets of correlations from standard 
texts have been prepared (in addition to the default correlations), allowing the user 
to *specify* which set to use.

## CCPS

The set of correlations for windspeed and atmospheric dispersion given in the CCPS 
*Guidelines for Consequence Analysis of Chemical Releases* (AIChE/CCPS 1999) and other
CCPS publications for quantitative risk assessment. There are two equation sets given 
by the CCPS: one for rural terrain and one for urban terrain.

+ `CCPSRural`
+ `CCPSUrban`

!!! note "Corrections"
    The correlations given for the vertical plume dispersion, $\sigma_z$, in urban 
    terrain in the CCPS references contains two typos. These have been corrected as per
    Griffiths (1994)


## ISC3

The set of correlations for windspeed and plume dispersion given in the ISC3 model 
description (EPA 1995). There are two equation sets: one for rural terrain and one for 
urban terrain.

+ `ISC3Rural`
+ `ISC3Urban`


## TNO

The set of correlations for atmospheric dispersion given in the TNO *Yellow Book* 
(Bakkum and Duijm 2005) for both continuous (plume) and instantaneous (puff) passive 
releases.

+ `TNO` 


## Turner

The standard plume dispersion plots from Turner (1970) as presented in Lees (1996)

+ `Turner`
