# GasDispersion.jl
[![LICENSE](https://img.shields.io/badge/license-MIT-lightgrey.svg)](https://github.com/aefarrell/GasDispersion.jl/blob/main/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-dev-blue)](https://aefarrell.github.io/GasDispersion.jl/dev/)
[![Build Status](https://github.com/aefarrell/GasDispersion.jl/workflows/CI/badge.svg)](https://github.com/aefarrell/GasDispersion.jl/actions)
[![Coverage](https://codecov.io/gh/aefarrell/GasDispersion.jl/branch/main/graph/badge.svg?token=PB3LOR80K2)](https://codecov.io/gh/aefarrell/GasDispersion.jl)

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

This scenario is adapted from CCPS *Guidelines for Consequence Analysis of
Chemical Releases*, CCPS, pg 47.

Suppose we wish to model the dispersion of gaseous propane from a leak from a storage tank, where the leak is from a 10 mm hole that is 3.5 m above the ground and the propane is at 25°C and 4barg. Assume the discharge coefficient $c_{D} = 0.85$

For ambient conditions we assume the atmosphere is dry air at standard conditions 
of 1atm and 25°C, with a windspeed of 1.5m/s and class F stability (a "worst case"
atmospheric stability), the default atmosphere if nothing else is specified.


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
       k = 1.15,             # heat capacity ratio, from Crane's
       temperature = 298.15, # K
       pressure = 501325,    # Pa
       height = 3.5)         # m, height of hole above the ground
```

This generates a `Scenario` defined for a gas jet discharging into dry air
at standard conditions. Once we have this defined we can determine the
concentration at any point downwind of the release point, assuming the release
is a continuous plume, using

```julia
# returns a callable
p = plume(scn, GaussianPlume)

p(x,y,z) # gives the concentration in kg/m^3 at the point x, y, z
```