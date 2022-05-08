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

Suppose a chemical release of some substance with a release rate of 1kg/s, at a
height of 1m. Using some standard engineering estimates we might end up with a
release scenario with the following parameters:

```julia
using GasDispersion

s = Scenario(
    1.0,   # mass emission rate, kg/s
    10.0,  # release duration, s
    0.25,  # jet diameter, m
    15.67, # jet velocity, m/s
    1.3,   # jet density, kg/m^3
    101325,# release_pressure, Pa
    450,   # release temperature, K
    1.0,   # release height, m
    1.5,   # windspeed, m/s
    1.225, # ambient density, kg/m^3
    101325,# ambient pressure, Pa
    298.15,# ambient temperature, K
    "F",   # pasquill stability class
)
```
Once we have this defined we can determine the concentration at any point
downwind of the release point, assuming the release is a continuous plume, using

```julia
# returns a function
c = plume(s)

c(x,y,z) # gives the concentration in kg/m^3 at the point x, y, z
```

where the coordinate system is such that the release point is at x=0, y=0, z=h

Similarly we could model an instantaneous release, assuming all of the mass was
released during 1 second, using a "puff" model
```julia
# returns a function
c = puff(s)

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

```@docs
Scenario
```

A `scenario_builder` function exists to help create valid `Scenario`s for
various standard release scenarios.
```@docs
scenario_builder
```


## Plume Models

Plume models are models of continuous, steady-state, releases and are time
independent, this includes, for example, emissions from elevated stacks.

```@docs
plume
```


## Puff Models

Puff models are for "instantaneous" releases or other time-dependent releases,
this often includes, for example, releases of vapour clouds.

```@docs
puff
```


## Future

In the future I would like to implement more dense gas dispersion models, as
well as a more diverse set of release scenarios.
