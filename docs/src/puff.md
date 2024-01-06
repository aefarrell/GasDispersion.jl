# Puff Models

Puff models are for "instantaneous" releases or other time-dependent releases.

```@docs
puff
```

## Gaussian Puffs

```@docs
puff(::Scenario, ::Type{GaussianPuff})
```

A gaussian puff model assumes the release is instantaneous, and all mass is concentrated in a single point. The cloud then disperses as it moves downwind with the concentration profile is given by a series of gaussians with dispersions $\sigma_x$, $\sigma_y$, and $\sigma_z$, which are found from correlations tabulated per stability class. Similarly to the plume model, a ground reflection term is included to correct for the fact that material cannot pass through the ground.

```math
c_{puff} = { V_i \over { (2 \pi)^{3/2} \sigma_x \sigma_y \sigma_z } } 
\exp \left( -\frac{1}{2} \left( {x - ut} \over \sigma_x \right)^2 \right) 
\exp \left( -\frac{1}{2} \left( {y} \over \sigma_y \right)^2 \right)  \left[ \exp \left( -\frac{1}{2} \left( {z - h} \over \sigma_z \right)^2 \right) 
+ \exp \left( -\frac{1}{2} \left( {z + h} \over \sigma_z \right)^2 \right)\right]
```

with
-  $c$ - concentration, volume fraction
-  $V_i$ - volume of released material (m^3)
- *u* - windspeed (m/s)
-  $\sigma_x$  - downwind dispersion (m)
-  $\sigma_y$  - crosswind dispersion (m)
-  $\sigma_z$  - vertical dispersion (m)
- *h* - release elevation (m)

The model assumes the initial release is a single point, with no dimensions. Unlike the plume model, this concentration is a function of time.

### Downwind dispersion correlations

The downwind dispersion, $\sigma_{x}$ is a function of downwind distance of the cloud center, $x_c$, as well as stability class

```math
\sigma_{x} = \delta {x_c}^{\beta}
```

Where $\delta$ and $\beta$ are identical to those tabulated for the crosswind dispersion.

### Crosswind dispersion correlations

The crosswind dispersion, $\sigma_{y}$ is a function of downwind distance of the cloud center, $x_c$, as well as stability class

```math
\sigma_{y} = \delta {x_c}^{\beta}
```

Where $\delta$, $\beta$, and $\gamma$ are tabulated based on stability class([AIChE/CCPS 1999](references.md)):

| Stability Class | $\delta$ | $\beta$ |
|:---------------:|:--------:|:-------:|
|        A        |  0.18    |   0.92  |
|        B        |  0.14    |   0.92  |
|        C        |  0.10    |   0.92  |
|        D        |  0.06    |   0.92  |
|        E        |  0.04    |   0.92  |
|        F        |  0.02    |   0.89  |


### Vertical dispersion correlations

The vertical dispersion, $\sigma_{z}$ is a function of downwind distance of the cloud center, $x_c$, as well as stability class

```math
 \sigma_{z} = \delta {x_c}^{\beta} 
```

Where $\delta$ and $\beta$ are tabulated based on stability class([AIChE/CCPS 1999](references.md)):

| Stability Class | $\delta$ | $\beta$ |
|:---------------:|:--------:|:-------:|
|        A        |  0.60    |  0.75   |
|        B        |  0.53    |  0.73   |
|        C        |  0.34    |  0.71   |
|        D        |  0.15    |  0.70   |
|        E        |  0.10    |  0.65   |
|        F        |  0.05    |  0.61   |

### Example

Suppose we wish to model the dispersion of gaseous propane from a leak from a storage tank, where the leak is from a 10mm hole that is 3.5m above the ground and the propane is at 25°C and 4barg. Assume the discharge coefficient $c_{D} = 0.85$. Assume the leak occurs for 10s. This scenario is adapted from CCPS *Guidelines for Consequence Analysis of Chemical Releases*([AIChE/CCPS 1999](references.md), 47)

First we define the scenario

```jldoctest propaneleak; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
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
       height = 3.5,     # m, height of hole above the ground
       duration = 10)    # s, duration of leak


# output

Substance: propane
    MW: 0.044096 kg/mol
    P_v: GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0) Pa
    ρ_g: 1.864931992847327 kg/m^3
    ρ_l: 526.13 kg/m^3
    T_ref: 288.15 K
    P_ref: 101325.0 Pa
    k: 1.142
    T_b: 231.02 K
    Δh_v: 425740 J/kg
    Cp_g: 1678 J/kg/K
    Cp_l: 2520 J/kg/K
HorizontalJet release:
    ṁ: 0.08991798763471508 kg/s
    Δt: 10.0 s
    d: 0.01 m
    u: 208.10961399327573 m/s
    h: 3.5 m
    P: 288765.2212333958 Pa
    T: 278.3846872082166 K
    f_l: 0.0
SimpleAtmosphere atmosphere:
    P: 101325.0 Pa
    T: 298.15 K
    u: 1.5 m/s
    h: 10.0 m
    rh: 0.0 %
    stability: ClassF

```

And then pass it to the `puff` function

```jldoctest propaneleak; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
g = puff(scn, GaussianPuff)

# output

GasDispersion.GaussianPuffSolution{Float64, ClassF, DefaultSet}(Scenario{Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, 10.0, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF)), :gaussian, 0.16344890861567063, 3.5, 1.150112899011524, ClassF, DefaultSet)

```

Where `g` is a callable which returns the concentration (in vol fraction) at any point. For example suppose we are interested in the concentration at some point 100m downwind of the release, along the centerline (y=0) and at a height of 2m, amd 86s after the start of the release.

```jldoctest propaneleak; output = true, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
x = 100 # m, the downwind distance
u = 1.5 # m/s, the windspeed
t = x/u # s, time when the cloud center arrives

g(100,0,2,86)

# output

0.0011119744194086872

```

```@setup propaneleak
using ..GasDispersion

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
       height = 3.5,     # m, height of hole above the ground
       duration = 10)    # s, duration of leak

g = puff(scn, GaussianPuff)

```

```@example propaneleak
using Plots

plot(g, 86; xlims=(90,110), ylims=(-10,10), aspect_ratio=:equal)

```

Using the `@gif` macro we can animate the arrival of the puff, visualizing how the cloud expands as it moves

```@example propaneleak
t = 86

@gif for t′ in range(0.85*t,1.15*t, length=50)

plot(g, t′; xlims=(85,115), ylims=(-10,10), clims=(0,1.5e-3), aspect_ratio=:equal)
    
end
```


## Integrated Gaussian Puffs

The integrated Gaussian puff model treats the release as a sum of $n$
equally spaced Gaussian puffs, starting at $t = 0$ to $t = \Delta t$. The
default behaviour is to take the limit $n \to \infty$.

```@docs
puff(::Scenario, ::Type{IntPuff})
```

## Britter-McQuaid Model

The Britter-McQuaid model is an empirical correlation for dense cloud
dispersion. The model generates an interpolation function for the average cloud
concentration and the cloud is rendered as a cylinder.


```@docs
puff(::Scenario, ::Type{BritterMcQuaidPuff})
```

## SLAB Horizontal Jet Model

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