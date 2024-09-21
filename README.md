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

```
pkg> add GasDispersion
```


## Example usage

Suppose we wish to model the dispersion of gaseous propane from a leak from a storage tank, where the leak is from a 10mm hole that is 3.5m above the ground and the propane is at 25°C and 4barg. Assume the discharge coefficient $c_{D} = 0.85$. This scenario is adapted from CCPS *Guidelines for Consequence Analysis of Chemical Releases*

First we define the scenario

```julia
using GasDispersion

propane = Substance(name="propane",
              molar_weight=0.044096,     # kg/mol
              liquid_density=526.13,     # kg/m³
              k=1.142,
              boiling_temp=231.02,       # K
              latent_heat=425740,        # J/kg
              gas_heat_capacity=1678,    # J/kg/K
              liquid_heat_capacity=2520) # J/kg/K

Patm = 101325 # Pa
P1 = 4e5 + Patm # Pa
T1 = 25 + 273.15 # K

scn = scenario_builder(propane, JetSource; 
       phase = :gas,
       diameter = 0.01,  # m
       dischargecoef = 0.85,
       temperature = T1, # K
       pressure = P1,    # Pa
       height = 3.5)     # m, height of hole above the ground
```

This generates a `Scenario` defined for a gas jet discharging into dry air
at standard conditions. 

An alternative to supplying all of the substance properties, which can get
rather onerous, is to use the library of substances and equations of state
provided in [Clapeyron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl)

```julia
using GasDispersion, Clapeyron

propane = Substance(PropaneRef()) # uses an empirical EoS for propane

propane = Substance(PR("propane"; idealmodel=ReidIdeal)) # uses the Peng-Robinson EoS
```

Once we have this defined we can determine the
concentration at any point downwind of the release point, assuming the release
is a continuous plume, using

```julia
p = plume(scn, GaussianPlume)
```

Where `p` is a callable which returns the concentration (in vol fraction) at any point. For example at 100m downwind and at a height of 2m

```julia
p(100,0,2)

# output

0.0002006455298894473

```

See [the documentation](https://aefarrell.github.io/GasDispersion.jl/dev/) for details on the additional plume and puff models and worked examples.