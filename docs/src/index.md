# GasDispersion.jl

GasDispersion.jl is a set of tools for atmospheric dispersion modeling of
gaseous releases, such as might occur during an emergency at a chemical plant
or more  routinely from a stack. This is intended to be the level of disperson
modeling used support consequence analysis or QRA such as is described in *Lee's
Loss Prevention in the Process Industries* or the CCPS *Guidelines for
Consequence Analysis of Chemical Releases*.

## Installation

GasDispersion.jl can be installed using Julia's built-in package manager. In a
Julia session, enter the package manager mode by hitting `]`, then run the
command

```
pkg> add https://github.com/aefarrell/GasDispersion.jl
```


## Example usage

!!! note "Worked examples"
    This is just a brief overview, detailed worked examples are available as [jupyter notebooks](https://nbviewer.org/github/aefarrell/GasDispersion.jl/tree/main/examples/)

This scenario is adapted from CCPS *Guidelines for Consequence Analysis of
Chemical Releases*, CCPS, pg 47.

Suppose we wish to model the dispersion of gaseous propane from a leak from a storage 
tank, where the leak is from a 10 mm hole that is 3.5 m above the ground and the 
propane is at 25°C and 4barg. Assume the discharge coefficient $c_{D} = 0.85$

For ambient conditions we assume the atmosphere is dry air at standard conditions
of 1atm and 25°C, with a windspeed of 1.5m/s and class F stability (a "worst case"
atmospheric stability)


```julia
using GasDispersion

propane = Substance(name = :propane,
                    gas_density = 9.7505, # kg/m^3, NIST Webbook
                    liquid_density = 526.13, # kg/m^3, NIST Webbook
                    reference_temp= 298.15, # K
                    reference_pressure= 501325, # Pa
                    boiling_temp = 231.04, # K, NIST Webbook
                    latent_heat = 16.25/44.0956, # J/kg, NIST Webbook
                    gas_heat_capacity = 1.6849,    # J/kg/K, NIST Webbook
                    liquid_heat_capacity = 2.2460) # J/kg/K, NIST Webbook

scn = scenario_builder(propane, JetSource;
       phase = :gas,
       diameter = 0.01, # m
       dischargecoef = 0.85,
       k = 1.15,         # heat capacity ratio, from Crane's
       temperature = T1, # K
       pressure = P1,    # Pa
       height = 3.5,     # m, height of hole above the ground
       duration = 1)     # s, duration of release
```

This generates a `Scenario` defined for a gas jet discharging into dry air
at standard conditions. Once we have this defined we can determine the
concentration at any point downwind of the release point, assuming the release
is a continuous plume, using

```julia
# returns a callable
pl = plume(scn, GaussianPlume)

pl(x,y,z) # gives the concentration in kg/m^3 at the point x, y, z
```
where the coordinate system is such that the release point is at x=0, y=0, z=h

Similarly we could model an instantaneous release, assuming all of the mass was
released during 1 second, using a "puff" model
```julia
# returns a function
pf = puff(scn, GaussianPuff)

p(x,y,z,t) # gives the concentration in kg/m^3 at the point x, y, z and time t
```



## Building Scenarios

A `Scenario` is a container for all of the information that a model may need to
produce a solution. The intention is for the `Scenario` to be re-usable, so that
the user may run the same scenario with multiple models without much difficulty.
Models also have specific parameters, those are handled in the model itself.

A `scenario_builder` function exists to help create valid `Scenario`s for
various standard release scenarios.
```@docs
scenario_builder

scenario_builder(::Substance, ::Type{JetSource}, ::Atmosphere)
```


## Plume Models

Plume models are models of continuous, steady-state, releases and are time
independent, this includes, for example, emissions from elevated stacks.

```@docs
plume
```

### Gaussian Plumes

A Gaussian plume is a steady state plume defined by a Gaussian distribution in
the y and z directions and an exponential decay in the x direction. The
dispersion parameters are correlated to the downwind distance.

```@docs
plume(::Scenario, ::Type{GaussianPlume})
```

### Simple Jet Plumes

The simple jet plume is a steady state turbulent jet defined by a Gaussian
distribution in the y and z directions. It is similar to the Gaussian plume,
however in this case the momentum forming the plume comes entirely from the jet
and not the ambient windspeed.

```@docs
plume(::Scenario, ::Type{SimpleJet})
```

### Britter-McQuaid Model

The Britter-McQuaid model is an empirical correlation for dense plume
dispersion. The model generates an interpolation function for the centerline
concentration at the downwind distance x.


```@docs
plume(::Scenario, ::Type{BritterMcQuaidPlume})
```


## Puff Models

Puff models are for "instantaneous" releases or other time-dependent releases.

```@docs
puff
```

### Gaussian Puffs

A Gaussian puff model is defined by gaussian distributions in the x, y, and z
directions and travels downwind at the ambient windspeed. The dispersion
parameters of the gaussians are correlated with the downwind distance of cloud
center.

```@docs
puff(::Scenario, ::Type{GaussianPuff})
```

### Integrated Gaussian Puffs

The integrated Gaussian puff model treats the release as a sum of $n$
equally spaced Gaussian puffs, starting at $t = 0$ to $t = \Delta t$. The
default behaviour is to take the limit $n \to \infty$.

```@docs
puff(::Scenario, ::Type{IntPuff})
```

### Britter-McQuaid Model

The Britter-McQuaid model is an empirical correlation for dense cloud
dispersion. The model generates an interpolation function for the average cloud
concentration and the cloud is rendered as a cylinder.


```@docs
puff(::Scenario, ::Type{BritterMcQuaidPuff})
```

### SLAB Horizontal Jet Model

The SLAB horizontal jet model is derived from the SLAB software package developed
by Donald L. Ermak at Lawrence Livermore National Laboratory. The model numerically
integrates a set of conservation equations for the given domain, automatically
transitioning from a steady-state plume model to a transient-state puff model as
the release terminates. The result is a set of cloud parameters that are interpolated
as a function of downwind distance and time to calculate the final concentration.

The SLAB model uses it's own built in models for atmospheric parameters, such as
windspeed and dispersion.

```@docs
puff(::Scenario, ::Type{SLAB})
```


## Equation Sets

The models above each depend upon several parameters which are themselves, often, 
correlations. For any given parameter there are several different correlations in the 
literature. To make this more transparent, sets of correlations from standard texts 
have been prepared (in addition to the default correlations), allowing the user to 
*specify* which set to use.

### CCPS

The set of correlations for windspeed and atmospheric dispersion given in the CCPS 
*Guidelines for Consequence Analysis of Chemical Releases* (AIChE/CCPS 1999) and other
CCPS publications for quantitative risk assessment. There are two equation sets given 
by the CCPS: one for rural terrain and one for urban terrain.

+ `CCPSRural`
+ `CCPSUrban`

!!! note Corrections
    The correlations given for the vertical plume dispersion, $\sigma_z$, in urban terrain in the CCPS references contains two typos. These have been corrected as per Griffiths (1994)


### ISC3

The set of correlations for windspeed and plume dispersion given in the ISC3 model 
description (EPA 1995). There are two equation sets: one for rural terrain and one for 
urban terrain.

+ `ISC3Rural`
+ `ISC3Urban`


### TNO

The set of correlations for atmospheric dispersion given in the TNO *Yellow Book* 
(Bakkum and Duijm 2005) for both continuous (plume) and instantaneous (puff) passive 
releases.

+ `TNO` 


### Turner

The standard plume dispersion plots from Turner (1970) as presented in Lees (1996)

+ `Turner`


# References

+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
+ Bakkum, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
+ Briggs, Gary A. 1969. *Plume Rise* Oak Ridge: U.S. Atomic Energy Commission
+ Briggs, Gary A. 1973. *Diffusion Estimation for Small Emissions. Preliminary Report*. United States. https://doi.org/10.2172/5118833
+ Britter, Rex E. and J. McQuaid. 1988. *Workbook on the Dispersion of Dense Gases. HSE Contract Research Report No. 17/1988*
+ Ermak, Donald L. 1990. *User's Manual for SLAB: An Atmospheric Dispersion Model For Denser-Than-Air Releases* Lawrence Livermore National Laboratory
+ EPA. 1995. *User's Guide for the Industrial Source Complex (ISC3) Dispersion Models, vol 2*. United States Environmental Protection Agency EPA-454/B-95-003b
+ Griffiths, R. F. 1994. "Errors in the use of the Briggs parameterization for atmospheric dispersion coefficients." *Atmospheric Environment* 28(17):2861-2865 https://doi.org/10.1016/1352-2310(94)90086-8
+ Lees, Frank P. 1996. *Loss Prevention in the Process Industries, 2nd ed*. Oxford: Butterworth-Heinemann
+ Long, V.D. 1963. "Estimation of the Extent of Hazard Areas Round a Vent." *Chem. Process Hazard*. II:6
+ Pasquill, Frank. 1974. *Atmospheric Diffusion, 2nd Ed.* New York: Halstead Press, New York
+ Spicer, Thomas O. and Jerry A. Havens. 1988. *Development of Vapor Dispersion Models for Non-Neutrally Buoyant Gas Mixtures--Analysis of TFI/NH3 Test Data*. United States.
+ Seinfeld, John H. 1986. *Atmospheric Chemistry and Physics of Air Pollution*. New York: John Wiley and Sons
+ Turner, D. Bruce. 1970. *Workbook of Atmospheric Dispersion Estimates*. United States Environmental Protection Agency.
