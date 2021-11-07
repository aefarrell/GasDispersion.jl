# GasDispersion.jl

[![Build Status](https://github.com/aefarrell/GasDispersion.jl/workflows/CI/badge.svg)](https://github.com/aefarrell/GasDispersion.jl/actions)
[![Coverage](https://codecov.io/gh/aefarrell/GasDispersion.jl/branch/main/graph/badge.svg?token=PB3LOR80K2)](https://codecov.io/gh/aefarrell/GasDispersion.jl)

GasDispersion.jl aims to bring together several models for dispersion modeling
of chemical releases with a consistent interface.

## Example usage

Suppose a chemical release of some substance with a release rate of 1kg/s, at a
height of 1m. Using some standard engineering estimates we might end up with a
release scenario with the following parameters:

```julia
scenario = Scenario(
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
plume_conc = plume(scenario, model=:gaussian)

plume_conc(x,y,z) # gives the concentration in kg/m^3 at the point x, y, z
```

where the coordinate system is such that the release point is at x=0, y=0, z=h

Similarly we could model an instantaneous release, assuming all of the mass was
released during 1 second, using a "puff" model
```julia
# returns a function
puff_conc = puff(scenario, model=:gaussian)

puff_conc(x,y,z,t) # gives the concentration in kg/m^3 at the point x, y, z and time t
```

## Future

Currently the only models implemented are simple gaussian plumes and puffs, in
the future I would like to implement some dense gas models as well.

Additionally constructing a scenario is relatively user unfriendly, and I would
like to implement some helper functions to generate `Scenario`s given some
minimal input.
