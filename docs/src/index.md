# GasDispersion.jl

GasDispersion.jl is a set of tools for atmospheric dispersion modeling of
gaseous releases, such as might occur during an emergency at a chemical plant
or more  routinely from a stack. This is intended to be the level of disperson
modeling used support consequence analysis or QRA such as is described in *Lee's
Loss Prevention in the Process Industries* or the CCPS *Guidelines for
Consequence Analysis of Chemical Releases*.

```@contents
Pages = ["scenarios.md", "plume.md", "puff.md", "equation_sets.md", "function_index.md"]
Depth = 2
```

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
                    latent_heat = 425740.0, # J/kg, 
                    gas_heat_capacity = 1678.0, # J/kg/K, 
                    liquid_heat_capacity = 2520.0) # J/kg/K, 

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

pl(x,y,z) # gives the concentration in vol fraction at the point x, y, z
```
where the coordinate system is such that the release point is at x=0, y=0, z=h

Similarly we could model an instantaneous release, assuming all of the mass was
released during 1 second, using a "puff" model
```julia
# returns a function
pf = puff(scn, GaussianPuff)

p(x,y,z,t) # gives the concentration in vol fraction at the point x, y, z and time t
```

