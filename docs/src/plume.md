# Plume Models

Plume models are models of continuous, steady-state, releases and are time
independent, this includes, for example, emissions from elevated stacks.

```@docs
plume
```

## Gaussian Plumes

```@docs
plume(::Scenario{Substance,VerticalJet,Atmosphere}, ::Type{GaussianPlume})
```

The gaussian plume model assumes that concentration profile in the crosswind (y) and vertical (z) directions follow gaussian distributions with dispersions $\sigma_y$ and $\sigma_z$, respectively. This model can be derived from an advection-diffusion equation, or simply taken as a given.

The basic gaussian would have the plume expand downward beyond the ground, to correct for this an additional term for *ground reflection* is added. This is equivalent to adding a mirror image source reflected across the x-z plane, and causes mass to accumulate along the ground instead of simply disappearing (as would happen in the naive case).

The concentration is then given by:

```math

 c = {Q_i \over 2 \pi u \sigma_{y} \sigma_{z} }
\exp \left[ -\frac{1}{2} \left( y \over \sigma_{y} \right)^2 \right]
\left\{ \exp \left[ -\frac{1}{2} \left( { z -h } \over \sigma_{z} \right)^2 \right] + \exp \left[ -\frac{1}{2} \left( { z + h } \over \sigma_{z} \right)^2 \right] \right\} 

```

with
-  $c$  - concentration, volume fraction
-  $Q_i$  - volumetric emission rate of the species, m^3/s
- *u* - windspeed, m/s
-  $\sigma_y$  - crosswind dispersion, m
-  $\sigma_z$  - vertical dispersion, m
- *h* - release elevation, m

Three important parameters are determined from correlations, which are functions of the atmospheric stability: the windspeed at the release point, the crosswind dispersion, and the vertical dispersion.

### Crosswind dispersion correlations

The crosswind dispersion, $\sigma_{y}$ is a function of downwind distance, $x$ as well as stability class

```math

 \sigma_{y} = \delta x^{\beta} 

```

Where $\delta$, $\beta$, and $\gamma$ are tabulated based on stability class ([Spicer and Havens 1988](references.md), 112):

| Stability Class | $\delta$ | $\beta$ |
|:---------------:|:--------:|:-------:|
|        A        |  0.423   |   0.9   |
|        B        |  0.313   |   0.9   |
|        C        |  0.210   |   0.9   |
|        D        |  0.136   |   0.9   |
|        E        |  0.102   |   0.9   |
|        F        |  0.0674  |   0.9   |


### Vertical dispersion correlations

The vertical dispersion, $\sigma_{z}$ is a function of downwind distance, $x$ as well as stability class

```math

 \sigma_{z} = \delta x^{\beta} \exp \left( \gamma \left( \ln x \right)^2 \right)

```

Where $\delta$ and $\beta$ are tabulated based on stability class ([Seinfeld 1986](references.md))

| Stability Class | $\delta$ | $\beta$ | $\gamma$ |
|:---------------:|:--------:|:-------:|:--------:|
|        A        |  107.7   | -1.7172 |  0.2770  |
|        B        |  0.1355  |  0.8752 |  0.0136  |
|        C        |  0.09623 |  0.9477 | -0.0020  |
|        D        |  0.04134 |  1.1737 | -0.0316  |
|        E        |  0.02275 |  1.3010 | -0.0450  |
|        F        |  0.01122 |  1.4024 | -0.0540  |

### Example

Suppose we wish to model the dispersion of gaseous propane from a leak from a storage tank, where the leak is from a 10mm hole that is 3.5m above the ground and the propane is at 25°C and 4barg. Assume the discharge coefficient $c_{D} = 0.85$. This scenario is adapted from CCPS *Guidelines for Consequence Analysis of Chemical Releases*([AIChE/CCPS 1999](references.md), 47)

First we define the scenario

```jldoctest gaussplume; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
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
    Δt: Inf s
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

And then pass it to the `plume` function

```jldoctest gaussplume; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
g = plume(scn, GaussianPlume)

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.NoPlumeRise, ClassF, DefaultSet}(Scenario{Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF)), :gaussian, 0.9999999999999998, 0.01634489086156706, 1.150112899011524, 3.5, GasDispersion.NoPlumeRise(), ClassF, DefaultSet)

```

Where `g` is a callable which returns the concentration (in vol fraction) at any point. For example at 100m downwind and at a height of 2m

```jldoctest gaussplume; output = true, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
g(100,0,2)

# output

0.0002006455298894473

```

Which is ~200ppm (vol). Beyond simply having a number, we may want a plan-view of the plume at a given height, say 2m.

```@setup gaussplume
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
       height = 3.5)     # m, height of hole above the ground

g = plume(scn, GaussianPlume)
```

```@example gaussplume
using Plots

plot(g, xlims=(0,100), ylims=(-10,10), height=2)
```

We might want instead to look at the concentration at the release height: 3.5 m, zoom in, and look at concentrations in the vicinity of the LEL.

```@example gaussplume
LEL = 0.021 # v/v, LEL from CAMEO Chemicals
UEL = 0.095 # v/v, UEL from CAMEO Chemicals

plot(g, xlims=(0,100), ylims=(-5,5), height=3.5, clims=(0,LEL))
```

These plan views stretch out the crosswind distance, but we can change the aspect ratio, to give a sense of how skinny the plume actually is

```@example gaussplume
plot(g, xlims=(0,50), ylims=(-10,10), height=3.5, clims=(0,LEL),
     aspect_ratio=:equal)
```

## Simple Jet Plumes

```@docs
plume(::Scenario, ::Type{SimpleJet})
```

Simple jet dispersion models are a useful tool for evaluating dispersion near the region where a jet release is occurring. They are based on a simplified model where the air is stationary and all of the momentum needed to mix the release is supplied by the jet. This is in some ways the opposite assumptions than are used in the Gaussian Plume model -- where the release is assumed to have negligible velocity and the momentum is entirely supplied by the wind.

```math
c = k_2 c_0 \left( d \over z \right) \sqrt{ \rho_j \over \rho_a } \exp \left( - \left( k_3 { y \over x } \right)^2 \right) \left[ \exp \left( - \left( k_3 { (z-h) \over x }\right)^2 \right) + \exp \left( - \left( k_3 { (z+h) \over x }\right)^2 \right) \right]
```

with
+ *c* - concentration, volume fraction
+  $k_2$ and $k_3$ - model parameters
+  $c_0$ - initial concentration, volume fraction
+ *d* - diameter of the jet, m
+  $\rho_j$ - initial density of the jet material, kg/m^3
+  $\rho_a$ - density of the ambient atmosphere, kg/m^3

### Model Parameters

The model parameters $k_2$ and $k_3$ are per [Long (1963)](references.md)

|       |     |
|------:|:----|
| $k_2$ | 6.0 |
| $k_3$ | 5.0 |

the initial concentration is calculated from the mass flowrate and volumetric flowrate

```math
 c_0 = { Q_i \over Q } = { \dot{m} \over \rho Q } = { \dot{m} \over { \rho \frac{\pi}{4} d^2 u } } 
```

### Example

Suppose we wish to model the dispersion of gaseous propane using the same scenario, `scn`, worked out above.

```jldoctest gaussplume; output = true, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
j = plume(scn, SimpleJet)
j(100,0,2)

# output

0.002485496609730624

```

```@example gaussplume
j = plume(scn, SimpleJet) # hide
plot(j, xlims=(0,100), ylims=(-10,10), height=2)
```


## Britter-McQuaid Model

```@docs
plume(::Scenario, ::Type{BritterMcQuaidPlume})
```

The Britter-McQuaid model is based on the *Workbook on the Dispersion of Dense Gases*([Britter and McQuaid 1988](references.md)) which uses a series of correlations relating maximum center-line concentrations to downwind distances based upon actual releases. The model of the plume is a series of rectangular slices with constant, average, concentration throughout -- giving a "top-hat" model. Points outside the defined plume are assumed to have zero concentration.

The original workbook gives the correlations as simply curves on a page, this model uses digitizations of those correlations([AIChE/CCPS 1999](references.md), 118)

### Example

This example is based on the Burro LNG dispersion results, in which LNG was released at ground-level, as given in *Guidelines for Consequence Analysis* ([AIChE/CCPS 1999](references.md), 122):
- release temperature: -162°C
- release rate: 0.23 m³/s (liquid)
- release duration: 174 s
- windspeed at 10m: 10.9 m/s
- LNG liquid density (at release conditions): 425.6 kg/m³
- LNG gas density (at release conditions): 1.76 kg/m³

and we assume the atmosphere was otherwise at ambient conditions of 298K and 1atm.

```jldoctest burro; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
using GasDispersion

lng = Substance(name = :LNG,
              molar_weight = 0.01604, # kg/mol, Methane
              gas_density = 1.76,
              liquid_density = 425.6,
              reference_temp=(273.15-162),
              reference_pressure=101325.0,
              boiling_temp = 111.6, # Methane, NIST Webbook
              latent_heat = 509880.0,  # J/kg, Methane
              gas_heat_capacity = 2240.0, # J/kg/K, Methane
              liquid_heat_capacity = 3349.0) # J/kg/K, Methane

ṁ = 0.23*lng.ρ_l # liquid spill rate times liquid density
Q = ṁ/lng.ρ_g    # gas volumetric flowrate

d = 1
A = (π/4)*(d)^2 # release area, assuming a diameter of 1m.
u = Q/A         # initial jet velocity

r = HorizontalJet( mass_rate = ṁ,
                  duration = 174,
                  diameter = d,
                  velocity = u,
                  height = 0.0,
                  pressure = 101325.0,
                  temperature = (273.15-162),
                  fraction_liquid = 0.0)

a = SimpleAtmosphere(windspeed=10.9, temperature=298, stability=ClassF)

scn = Scenario(lng,r,a)

# output

Substance: LNG 
    MW: 0.01604 kg/mol 
    P_v: GasDispersion.Antoine{Float64}(8.814018574933064, 983.6444729625298, 0.0) Pa 
    ρ_g: 1.76 kg/m^3 
    ρ_l: 425.6 kg/m^3 
    T_ref: 111.14999999999998 K 
    P_ref: 101325.0 Pa 
    k: 1.4  
    T_b: 111.6 K 
    Δh_v: 509880.0 J/kg 
    Cp_g: 2240.0 J/kg/K 
    Cp_l: 3349.0 J/kg/K 
HorizontalJet release:
    ṁ: 97.888 kg/s 
    Δt: 174.0 s 
    d: 1.0 m 
    u: 70.81526849717933 m/s 
    h: 0.0 m 
    P: 101325.0 Pa 
    T: 111.14999999999998 K 
    f_l: 0.0  
SimpleAtmosphere atmosphere:
    P: 101325.0 Pa 
    T: 298.0 K 
    u: 10.9 m/s 
    h: 10.0 m 
    rh: 0.0 % 
    stability: ClassF  


```

Generating a solution using the Britter-McQuaid model is quite simple

```jldoctest burro; output = true, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
bm = plume(scn, BritterMcQuaidPlume)

bm(367,0,0)

# output

0.050717667650511944

```

From the reference, we expect the concentration to be ~5%, which is within acceptable margins given the imprecision that comes with using correlations (especially with pencil and paper).