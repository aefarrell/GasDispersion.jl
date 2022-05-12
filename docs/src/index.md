# GasDispersion.jl

GasDispersion.jl aims to bring together several models for dispersion modeling
of chemical releases with a consistent interface. Currently it is very much
under development and significant portions of the code and interface are subject
to change.

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
Chemical Releases*, CCPS, pg 40.

Suppose a leak of a liquid from a storage tank through a hole. At the hole, the
liquid pressure is 120kPa, the liquid has a density of 490 kg/m³, and is in
thermal equilibrium with the environment. Assume a circular hole with a
discharge coefficient of 0.63 and a diameter of 1cm.

For ambient conditions we assume the atmosphere is at standard conditions of
1atm and 25°C, with a windspeed of 1.5m/s and class F stability (a "worst case"
atmospheric stability)

Assumptions:
+ discharge coefficient, 0.63
+ diameter of the hole, 0.010m
+ liquid density, 490kg/m3
+ liquid temperature, 25C or 298.15K
+ liquid pressure, 120kPA or 120000Pa
+ ambient pressure, 1atm or 101325Pa
+ ambient temperature, 25C or 298.15K
+ density of air at standard conditions, 1.225kg/m3
+ ambient windspeed, 1.5m/s
+ Pasquill-Gifford stability F

```julia
using GasDispersion

source=JetSource(phase=:liquid, dischargecoef=0.63, diameter=0.01,
                 pressure=120000, temperature=298.15, density=490, height=1)

s = scenario_builder(source, Ambient())
```

Returns a `Scenario` defined for a liquid jet discharging into the air at
standard conditions. Once we have this defined we can determine the
concentration at any point downwind of the release point, assuming the release
is a continuous plume, using

```julia
# returns a function
c = plume(s, GaussianPlume())

c(x,y,z) # gives the concentration in kg/m^3 at the point x, y, z
```

where the coordinate system is such that the release point is at x=0, y=0, z=h

Similarly we could model an instantaneous release, assuming all of the mass was
released during 1 second, using a "puff" model
```julia
# returns a function
c = puff(s, GaussianPuff())

c(x,y,z,t) # gives the concentration in kg/m^3 at the point x, y, z and time t
```



## Building Scenarios

A `Scenario` is a container for all of the information that a model may need to
produce a solution. The intention is for the `Scenario` to be re-usable, so that
the user may run the same scenario with multiple models without much difficulty.
Models also have specific parameters, those are handled in the model itself.

By default a `Scenario` can have any field `missing`, this is because not all
models require all fields. Each model then verifies that none of the necessary
fields are `missing` and throws an error otherwise.

A `scenario_builder` function exists to help create valid `Scenario`s for
various standard release scenarios.
```@docs
scenario_builder

scenario_builder(::JetSource, ::Atmosphere)
```


## Plume Models

Plume models are models of continuous, steady-state, releases and are time
independent, this includes, for example, emissions from elevated stacks.

```@docs
plume
```

### Gaussian Plumes

A gaussian plume is a steady state plume defined by a gaussian distribution in
the y and z directions and an exponential decay in the x direction. The
dispersion parameters are correlated to the downwind distance.

```@docs
plume(::Scenario, ::GaussianPlume)
```

### Simple Jet Plumes

The simple jet plume is a steady state turbulent jet defined by a gaussian
distribution in the y and z directions. It is similar to the gaussian plume,
however in this case the momentum forming the plume comes entirely from the jet
and not the ambient windspeed.

```@docs
plume(::Scenario, ::SimpleJet)
```

### Britter-McQuaid Model

The Britter-McQuaid model is an empirical correlation for dense plume
dispersion. The model generates an interpolation function for the centerline
concentration at the downwind distance z.


```@docs
plume(::Scenario, ::BritterMcQuaidPlume)
```


## Puff Models

Puff models are for "instantaneous" releases or other time-dependent releases,
this often includes, for example, releases of vapour clouds.

```@docs
puff
```

### Gaussian Puffs

A gaussian puff model is defined by gaussian distributions in the x, y, and z
directions and travels downwind at the ambient windspeed. The dispersion
parameters of the gaussians are correlated with the downwind distance.

```@docs
puff(::Scenario, ::GaussianPuff)
```


## Future

In the future I would like to implement more dense gas dispersion models, as
well as a more diverse set of release scenarios.
