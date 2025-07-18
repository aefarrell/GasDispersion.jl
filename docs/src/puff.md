# Puff Models

Puff models are for "instantaneous" releases or other time-dependent releases.

```@docs
puff
```

## Gaussian Puff Models

### Simple Gaussial Puffs

```@docs
puff(::Scenario, ::GaussianPuff)
```

A simple gaussian puff model assumes the release is instantaneous, and all mass is concentrated in a single point. The cloud then disperses as it moves downwind with the concentration profile is given by a series of gaussians with dispersions $\sigma_x$, $\sigma_y$, and $\sigma_z$, which are found from correlations tabulated per stability class. Similarly to the plume model, a ground reflection term is included to correct for the fact that material cannot pass through the ground.

```math
c_{puff} = { m_i \over { (2 \pi)^{3/2} \sigma_x \sigma_y \sigma_z } } 
\exp \left( -\frac{1}{2} \left( {x - ut} \over \sigma_x \right)^2 \right) 
\exp \left( -\frac{1}{2} \left( {y} \over \sigma_y \right)^2 \right)  \\
\left[ \exp \left( -\frac{1}{2} \left( {z - h} \over \sigma_z \right)^2 \right) 
+ \exp \left( -\frac{1}{2} \left( {z + h} \over \sigma_z \right)^2 \right)\right]
```

with
- *c* - concentration, kg/m^3
-  $m_i$ - mass of released material, kg
- *u* - windspeed, m/s
-  $\sigma_x$  - downwind dispersion, m
-  $\sigma_y$  - crosswind dispersion, m
-  $\sigma_z$  - vertical dispersion, m
- *h* - release elevation, m

The model assumes the initial release is a single point, with no dimensions. Unlike the plume model, this concentration is a function of time. The model converts the final concentration to volume fraction, assuming the puff is a gas at ambient conditions.

#### Downwind dispersion correlations

The downwind dispersion, $\sigma_{x}$ is a function of downwind distance of the cloud center, $x_c$, as well as stability class

```math
\sigma_{x} = \delta {x_c}^{\beta}
```

Where $\delta$ and $\beta$ are identical to those tabulated for the crosswind dispersion.

#### Crosswind dispersion correlations

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


#### Vertical dispersion correlations

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

#### Example

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

scn = scenario_builder(propane, JetSource(); 
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
    stability: ClassF()

```

And then pass it to the `puff` function

```jldoctest propaneleak; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
g = puff(scn, GaussianPuff())

# output

GasDispersion.GaussianPuffSolution{Float64, BasicEquationSet{DefaultWind, CCPSPuffσx, CCPSPuffσy, CCPSPuffσz}}(Scenario{Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, 10.0, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF())), :gaussian, 0.8991798763471508, 1.8023818673116125, 3.5, 1.150112899011524, BasicEquationSet{DefaultWind, CCPSPuffσx, CCPSPuffσy, CCPSPuffσz}(DefaultWind(), CCPSPuffσx(), CCPSPuffσy(), CCPSPuffσz()))

```

Where `g` is a callable which returns the concentration (in vol fraction) at any point. For example suppose we are interested in the concentration at some point 100m downwind of the release, along the centerline (y=0) and at a height of 2m, amd 86s after the start of the release.

```jldoctest propaneleak; output = true, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
g(100,0,2,86)

# output

0.003394005492341503

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

scn = scenario_builder(propane, JetSource(); 
       phase = :gas,
       diameter = 0.01,  # m
       dischargecoef = 0.85,
       temperature = T1, # K
       pressure = P1,    # Pa
       height = 3.5,     # m, height of hole above the ground
       duration = 10)    # s, duration of leak

g = puff(scn, GaussianPuff())

```

```@example propaneleak
using Plots

plot(g, 86; xlims=(90,110), ylims=(-10,10), aspect_ratio=:equal)

```

Using the `@gif` macro we can animate the arrival of the puff, visualizing how the cloud expands as it moves

```@example propaneleak
t = 86

@gif for t′ in range(0.85*t,1.15*t, length=50)

plot(g, t′; xlims=(85,115), ylims=(-10,10), clims=(0,5e-3), aspect_ratio=:equal)
    
end
```

### Palazzi Short Duration Puff Model
```@docs
    puff(::Scenario, ::Palazzi)
```
The `Palazzi` model integrates over the Gaussian puff model for a short duration release ([Palazzi 1982](references.md)), where the dispersion parameters $\sigma$ are assumed to be independent time.

```math
c\left(x,y,z,t\right) = \chi\left(x,y,z\right) \frac{1}{2} \left( \mathrm{erf} \left( { {x - u (t-\Delta t)} \over \sqrt{2} \sigma_x } \right) - \mathrm{erf} \left( { {x - u t} \over \sqrt{2} \sigma_x } \right)  \right)
```

where $\chi$ is a Gaussian plume model and $\sigma_x$ is the downwind dispersion parameter. The model assumes the initial release is a single point, with no dimensions. Additionally, the model converts the final concentration to volume fraction, assuming the puff is a gas at ambient conditions.

By default the Palazzi model assumes a simple Gaussian plume model, $\chi$, for which this is a "correction", and uses plume dispersion parameters with $\sigma_x \left( x \right) = \sigma_y \left( x \right)$. However the user is free to use *any* plume model which is a subtype of `Plume`, and any equation set that implements downwind dispersion.

There are multiple variations of the Palazzi short duration model, depending on how the downwind dispersion, $\sigma_x$ , is calculated:
- `:default` follows Palazzi and calculates $\sigma_x$  at the downwind distance *x*
- `:intpuff` calculates $\sigma_x$  at the downwind distance to the cloud center at the start and end of the cloud, $ut$ and $u \left(t-\Delta t\right)$
- `:tno` follows the TNO Yellow Book eqn 4.60b, using the distance *x* while the plume is still attached to the release point, and the distance to the cloud center, *ut*, afterwards


### Integrated Gaussian Puff Model

```@docs
puff(::Scenario, ::IntPuff)
```

The `IntPuff` model treats a release as a sequence of $n$ gaussian puffs, each one corresponding to $\frac{1}{n}$ of the total mass of the release.

```math
c\left(x,y,z,t\right) = \sum_{i}^{n-1} { {m_i \Delta t} \over n } { { \exp \left( -\frac{1}{2} \left( {x - u \left( t - i \delta t \right) } \over \sigma_x \right)^2 \right) } \over { \sqrt{2\pi} \sigma_x } } { { \exp \left( -\frac{1}{2} \left( {y} \over \sigma_y \right)^2 \right) } \over { \sqrt{2\pi} \sigma_y } }\\ \times { { \exp \left( -\frac{1}{2} \left( {z - h} \over \sigma_z \right)^2 \right) + \exp \left( -\frac{1}{2} \left( {z + h} \over \sigma_z \right)^2 \right) } \over { \sqrt{2\pi} \sigma_z } }
```


with 
- *c* - concentration, kg/m^3
-  $m_i$ - emission rate, kg/s
-  $\Delta t$ - total duration, s
-  $ \delta t = {\Delta t \over n} $ - puff interval, s
- *u* - windspeed, m/s
-  $\sigma_x$ - downwind dispersion, m
-  $\sigma_y$ - crosswind dispersion, m
-  $\sigma_z$ - downwind dispersion, m

The model assumes the initial release is a single point, with no dimensions. Additionally, model converts the final concentration to volume fraction, assuming the puff is a gas at ambient conditions. If no number, *n*, is specified the model defaults to an integral approximation similar to the [Palazzi Short Duration Puff Model](@ref).

#### Example

Continuing with the propane leak example from above, we now model the release as a sequence of 100 gaussian puffs. Essentially chopping the 10s over which the release happens into 0.1s intervals and releasing one puff per interval at a time for 10s.

```jldoctest propaneleak; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
ig = puff(scn, IntPuff(); n=100)

# output

GasDispersion.IntPuffSolution{Float64, Int64, BasicEquationSet{DefaultWind, CCPSPuffσx, CCPSPuffσy, CCPSPuffσz}}(Scenario{Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, 10.0, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF())), :intpuff, 0.08991798763471508, 1.8023818673116125, 10.0, 3.5, 1.150112899011524, 100, BasicEquationSet{DefaultWind, CCPSPuffσx, CCPSPuffσy, CCPSPuffσz}(DefaultWind(), CCPSPuffσx(), CCPSPuffσy(), CCPSPuffσz()))
```

At the same point as above the concentration has dropped
```jldoctest propaneleak; output = true, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
ig(100,0,2,86)

# output

0.0002521339225936648

```

```@setup propaneleak
ig = puff(scn, IntPuff; n=100)
```

```@example propaneleak

plot(ig, 86; xlims=(90,110), ylims=(-10,10), aspect_ratio=:equal)

```

For short duration releases the model approximates the integral when $n \to \infty$ this is the default behaviour or when `n=Inf`

```@example propaneleak
ig_inf = puff(scn, IntPuff)

plot(ig_inf, 86; xlims=(90,110), ylims=(-10,10), aspect_ratio=:equal)
```

## Box Models
### Britter-McQuaid Model

```@docs
puff(::Scenario, ::BritterMcQuaidPuff)
```

The Britter-McQuaid model is an empirical correlation for dense cloud
dispersion. The model generates an interpolation function for the average cloud
concentration and the cloud is rendered as a cylinder. The only correlations used
in the provided equationset are for windspeed.


## Integral Models
### SLAB Jet Model

```@docs
puff(::Scenario{Substance,VerticalJet,Atmosphere}, ::SLAB)
```

The SLAB jet model is derived from the SLAB software package developed by 
Donald L. Ermak at Lawrence Livermore National Laboratory. The model numerically
integrates a set of conservation equations for the given domain, automatically
transitioning from a steady-state plume model to a transient-state puff model as
the release terminates. The result is a set of cloud parameters that are interpolated
as a function of downwind distance and time to calculate the final concentration.

The SLAB model uses it's own built in models for atmospheric parameters, such as
windspeed and dispersion.