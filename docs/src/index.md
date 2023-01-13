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

```julia
pkg> add https://github.com/aefarrell/GasDispersion.jl
```


## Example usage

!!! note "Worked examples"
    This is just a brief overview, detailed worked examples are available as [jupyter notebooks](https://nbviewer.org/github/aefarrell/GasDispersion.jl/tree/main/examples/)

This scenario is adapted from CCPS *Guidelines for Consequence Analysis of
Chemical Releases*, CCPS, pg 47.

Suppose we wish to model the dispersion of gaseous propane from a leak from a storage tank, where the leak is from a 10 mm hole that is 3.5 m above the ground and the propane is at 25°C and 4barg. Assume the discharge coefficient $c_{D} = 0.85$

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
