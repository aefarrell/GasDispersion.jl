# GasDispersion.jl
[![LICENSE](https://img.shields.io/badge/license-MIT-lightgrey.svg)](https://github.com/aefarrell/GasDispersion.jl/blob/main/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-dev-blue)](https://aefarrell.github.io/GasDispersion.jl/dev/)
[![Build Status](https://github.com/aefarrell/GasDispersion.jl/workflows/CI/badge.svg)](https://github.com/aefarrell/GasDispersion.jl/actions)
[![Coverage](https://codecov.io/gh/aefarrell/GasDispersion.jl/branch/main/graph/badge.svg?token=PB3LOR80K2)](https://codecov.io/gh/aefarrell/GasDispersion.jl)

GasDispersion.jl aims to bring together several models for dispersion modeling
of chemical releases with a consistent interface.

## Installation

GasDispersion.jl can be installed using Julia's built-in package manager. In a
Julia session, enter the package manager mode by hitting `]`, then run the
command

```julia
pkg> add https://github.com/aefarrell/GasDispersion.jl
```


## Example usage

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

s = scenario_builder(120000, 298.15; windspeed=1.5, stability="F",
                     model=:jet, phase=:liquid, liquid_density=490,
                     hole_diameter=0.01, discharge_coeff=0.63)
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

## Future

In the future I would like to implement more dense gas dispersion models, as
well as a more diverse set of release scenarios.
