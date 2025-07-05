# Equation Sets

The various dispersion models each depend upon several parameters which are themselves, often, correlations. For any given parameter there are several different correlations in the literature. To make this more transparent, sets of correlations from standard texts have been prepared (in addition to the default correlations), allowing the user to *specify* which set to use.

There are seven equation sets for plume models which define the correlations for windspeed, crosswind dispersion, and vertical dispersion:
+ `DefaultSet` -- see [the defition in plume models](plume.md)
+ `CCPSRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSUrban` -- [AIChE/CCPS 1999](references.md)
+ `ISC3Rural` -- [US EPA 1995](references.md)
+ `ISC3Urban` -- [US EPA 1995](references.md)
+ `TNOPlume` -- *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))
+ `Turner` -- Dispersion parameters from [Turner (1970)](references.md) and [Lees 1996](references.md), windspeed is the default

Mostly these reference a smaller set of power-law wind correlations and dispersion correlations. The specifc details are given below with more details on the particular correlations in the corresponding sections below.

| Equation Set | Wind           | $\sigma_x$ | $\sigma_y$     | $\sigma_z$     | 
|:------------:|:--------------:|:----------:|:--------------:|:--------------:|
| `DefaultSet` | `DefaultWind`  | `Nothing`  | `Defaultσy`    | `Defaultσz`    |
| `CCPSRural`  | `IrwinRural`   | `Nothing`  | `BriggsRuralσy`| `BriggsRuralσz`|
| `CCPSUrban`  | `IrwinUrban`   | `Nothing`  | `BriggsUrbanσy`| `BriggsUrbanσz`|
| `ISC3Urban`  | `ISC3UrbanWind`| `Nothing`  | `BriggsUrbanσy`| `BriggsUrbanσz`|
| `ISC3Rural`  | `IrwinRural`   | `Nothing`  | `ISC3Ruralσy`  | `ISC3Ruralσz`  |
| `TNOPlume`   | `TNOWind`      | `Nothing`  | `TNOPlumeσy`   | `TNOPlumeσz`   |
| `Turner`     | `DefaultWind`  | `Nothing`  | `Turnerσy`     | `Turnerσz`     |

There are four equation sets for puff models which define the correlations for windspeed, downwind dispersion, crosswind dispersion, and vertical dispersion:
+ `DefaultPuffSet` -- see [the defition in plume models](plume.md)
+ `CCPSPuffRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSPuffUrban` -- [AIChE/CCPS 1999](references.md)
+ `TNOPuff` -- *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))

| Equation Set     | Wind           | $\sigma_x$    | $\sigma_y$  | $\sigma_z$     | 
|:----------------:|:--------------:|:-------------:|:-----------:|:--------------:|
| `DefaultPuffSet` | `DefaultWind`  | `CCPSPuffσx`  | `CCPSPuffσy`| `CCPSPuffσz`   |
| `CCPSPuffRural`  | `IrwinRural`   | `CCPSPuffσx`  | `CCPSPuffσy`| `CCPSPuffσz`   |
| `CCPSPuffUrban`  | `IrwinUrban`   | `CCPSPuffσx`  | `CCPSPuffσy`| `CCPSPuffσz`   |
| `TNOPuff`        | `TNOWind`      | `TNOPuffσz`   | `TNOPuffσy` | `TNOPuffσz`    |

An equation set intended for puff models can be used for a plume model, but not vice-versa unless otherwise noted. Equation sets for plume models do not typically define a downwind dispersion, $\sigma_x$, and without some additional details on how to handle that there is no way to calculate puff dispersion.

It is possible to define one's own equation set by either mixing and matching existing correlations. For example, suppose I want to use the TNO puff dispersion correlations but with the Irwin rural powerlaw wind profile:

```@example
using GasDispersion

MyPuffSet = BasicEquationSet(IrwinRural(),TNOPuffσz(),TNOPuffσy(),TNOPuffσz())
```

All of the pre-defined equation sets are simply constants set to a particular instance of `BasicEquationSet`. Additional correlations can be added for windspeed or dispersion, by overloading the internal functions `windspeed`, `downwind_dispersion`, `crosswind_dispersion`, and `vertical_dispersion`, but this can be dangerous as the internals are subject to change without notice. A better choice is to open an issue or pull request through GitHub, and it can be added to the next release.

### Example Usage

Using the same example scenario as the basic gaussian plume, we can explore the sensitivity to choice of model parameters using equation sets. Starting with the scenario definition:

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
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
    stability: ClassF()

```

The plume using the default equation set is simply this

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
dflt = plume(scn, GaussianPlume(), DefaultSet)

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.SimpleCrossTerm, GasDispersion.SimpleVerticalTerm, GasDispersion.NoPlumeRise, BasicEquationSet{DefaultWind, Nothing, Defaultσy, Defaultσz}, GasDispersion.ProblemDomain{Float64}}(Scenario{Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF())), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 1.150112899011524, 3.5, GasDispersion.SimpleCrossTerm(), GasDispersion.SimpleVerticalTerm(), GasDispersion.NoPlumeRise(), BasicEquationSet{DefaultWind, Nothing, Defaultσy, Defaultσz}(DefaultWind(), nothing, Defaultσy(), Defaultσz()), GasDispersion.ProblemDomain{Float64}(0.0, Inf, -Inf, Inf, 0.0, Inf))

```

For each of the plume equation sets we can easily create corresponding plume solutions

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
ccps_rurl = plume(scn, GaussianPlume(), CCPSRural)

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.SimpleCrossTerm, GasDispersion.SimpleVerticalTerm, GasDispersion.NoPlumeRise, BasicEquationSet{IrwinRural, Nothing, BriggsRuralσy, BriggsRuralσz}, GasDispersion.ProblemDomain{Float64}}(Scenario{Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF())), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 0.8420321686971456, 3.5, GasDispersion.SimpleCrossTerm(), GasDispersion.SimpleVerticalTerm(), GasDispersion.NoPlumeRise(), BasicEquationSet{IrwinRural, Nothing, BriggsRuralσy, BriggsRuralσz}(IrwinRural(), nothing, BriggsRuralσy(), BriggsRuralσz()), GasDispersion.ProblemDomain{Float64}(0.0, Inf, -Inf, Inf, 0.0, Inf))

```

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
ccps_urb = plume(scn, GaussianPlume(), CCPSUrban)

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.SimpleCrossTerm, GasDispersion.SimpleVerticalTerm, GasDispersion.NoPlumeRise, BasicEquationSet{IrwinUrban, Nothing, BriggsUrbanσy, BriggsUrbanσz}, GasDispersion.ProblemDomain{Float64}}(Scenario{Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF())), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 0.7989729675905327, 3.5, GasDispersion.SimpleCrossTerm(), GasDispersion.SimpleVerticalTerm(), GasDispersion.NoPlumeRise(), BasicEquationSet{IrwinUrban, Nothing, BriggsUrbanσy, BriggsUrbanσz}(IrwinUrban(), nothing, BriggsUrbanσy(), BriggsUrbanσz()), GasDispersion.ProblemDomain{Float64}(0.0, Inf, -Inf, Inf, 0.0, Inf))

```

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
isc3_rurl = plume(scn, GaussianPlume(), ISC3Rural)

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.SimpleCrossTerm, GasDispersion.SimpleVerticalTerm, GasDispersion.NoPlumeRise, BasicEquationSet{IrwinRural, Nothing, ISC3Ruralσy, ISC3Ruralσz}, GasDispersion.ProblemDomain{Float64}}(Scenario{Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF())), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 0.8420321686971456, 3.5, GasDispersion.SimpleCrossTerm(), GasDispersion.SimpleVerticalTerm(), GasDispersion.NoPlumeRise(), BasicEquationSet{IrwinRural, Nothing, ISC3Ruralσy, ISC3Ruralσz}(IrwinRural(), nothing, ISC3Ruralσy(), ISC3Ruralσz()), GasDispersion.ProblemDomain{Float64}(0.0, Inf, -Inf, Inf, 0.0, Inf))

```

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
isc3_urb = plume(scn, GaussianPlume(), ISC3Urban)

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.SimpleCrossTerm, GasDispersion.SimpleVerticalTerm, GasDispersion.NoPlumeRise, BasicEquationSet{ISC3UrbanWind, Nothing, BriggsUrbanσy, BriggsUrbanσz}, GasDispersion.ProblemDomain{Float64}}(Scenario{Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF())), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 1.0947417281650496, 3.5, GasDispersion.SimpleCrossTerm(), GasDispersion.SimpleVerticalTerm(), GasDispersion.NoPlumeRise(), BasicEquationSet{ISC3UrbanWind, Nothing, BriggsUrbanσy, BriggsUrbanσz}(ISC3UrbanWind(), nothing, BriggsUrbanσy(), BriggsUrbanσz()), GasDispersion.ProblemDomain{Float64}(0.0, Inf, -Inf, Inf, 0.0, Inf))

```

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
tno = plume(scn, GaussianPlume(), TNOPlume)

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.SimpleCrossTerm, GasDispersion.SimpleVerticalTerm, GasDispersion.NoPlumeRise, BasicEquationSet{TNOWind, Nothing, TNOPlumeσy, TNOPlumeσz}, GasDispersion.ProblemDomain{Float64}}(Scenario{Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF())), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 0.8753751236458281, 3.5, GasDispersion.SimpleCrossTerm(), GasDispersion.SimpleVerticalTerm(), GasDispersion.NoPlumeRise(), BasicEquationSet{TNOWind, Nothing, TNOPlumeσy, TNOPlumeσz}(TNOWind(), nothing, TNOPlumeσy(), TNOPlumeσz()), GasDispersion.ProblemDomain{Float64}(0.0, Inf, -Inf, Inf, 0.0, Inf))

```

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
turner = plume(scn, GaussianPlume(), Turner)

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.SimpleCrossTerm, GasDispersion.SimpleVerticalTerm, GasDispersion.NoPlumeRise, BasicEquationSet{DefaultWind, Nothing, Turnerσy, Turnerσz}, GasDispersion.ProblemDomain{Float64}}(Scenario{Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, Float64, GasDispersion.Antoine{Float64}, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF())), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 1.150112899011524, 3.5, GasDispersion.SimpleCrossTerm(), GasDispersion.SimpleVerticalTerm(), GasDispersion.NoPlumeRise(), BasicEquationSet{DefaultWind, Nothing, Turnerσy, Turnerσz}(DefaultWind(), nothing, Turnerσy(), Turnerσz()), GasDispersion.ProblemDomain{Float64}(0.0, Inf, -Inf, Inf, 0.0, Inf))

```

All of these plumes can then be plotted, to better visualize what is going on. These are identical plume models with the only differences being the windspeed correlation and the dispersion correlations.

```@eval
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
       height = 3.5)     # m, height of hole above the ground

dflt = plume(scn, GaussianPlume(), DefaultSet)
ccps_rurl = plume(scn, GaussianPlume(), CCPSRural)
ccps_urb = plume(scn, GaussianPlume(), CCPSUrban)
isc3_rurl = plume(scn, GaussianPlume(), ISC3Rural)
isc3_urb = plume(scn, GaussianPlume(), ISC3Urban)
tno = plume(scn, GaussianPlume(), TNOPlume)
turner = plume(scn, GaussianPlume(), Turner)

solset = [("Default", dflt), ("CCPSRural", ccps_rurl), 
          ("CCPSUrban", ccps_urb), ("ISC3Rural", isc3_rurl),
          ("ISC3Urban", isc3_urb), ("TNO", tno), ("Turner", turner)]

using Plots

xs = 1:1:500 
ys = range(-50,50,200)

plt1 = plot()    
plot!(plt1, xlabel="Downwind distance (m)", ylabel="Concentration C₃H₈ (ppmv)", title="Downwind concentration at groundlevel", legend=:topright)

for sol in solset
    lbl, pl = sol
    plot!(plt1, xs, pl.(xs,0,0).*1e6, label=lbl)
end

savefig(plt1,"example_fig-1.svg")

nothing
```

![](example_fig-1.svg)

## Windspeed

The most common windspeed profile is a power-law relationship:
```math
 u = u_{R} \left( z \over z_{R} \right)^{p}
```

There are four power-law correlations for windspeed:
+ `DefaultWind` -- see [the defition in release scenarios](scenarios.md)
+ `IrwinRural` -- [Irwin 1979](references.md)
+ `IrwinUrban` -- [Irwin 1979](references.md)
+ `ISC3UrbanWind` -- [US EPA 1995](references.md)

| Stability Class | `DefaultWind`  | `IrwinRural`  | `IrwinUrban`  | `ISC3UrbanWind` |
|:---------------:|:--------------:|:-------------:|:-------------:|:---------------:|
|        A        |  0.108         |  0.07         |  0.15         | 0.15            |
|        B        |  0.112         |  0.07         |  0.15         | 0.15            |
|        C        |  0.120         |  0.10         |  0.20         | 0.20            |
|        D        |  0.142         |  0.15         |  0.25         | 0.25            |
|        E        |  0.203         |  0.35         |  0.40         | 0.30            |
|        F        |  0.253         |  0.55         |  0.60         | 0.30            |

!!! note "Note"
    The ISC3Urban correlation is the same as the IrwinUrban except for stable atmospheres (class E and F) 

There are two correlations which uses a logarithmic profile based on Monin-Obukhov similarity theory.
+ `TNOWind` -- *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md)), with a default surface roughnes of 0.1m
+ `BusingerWind` -- [Businger *et al.* 1971](references.md)

!!! note "Surface Roughness"
    The `SimpleAtmosphere` type does not define a surface roughness, however wind profiles based on Monin-Obukhov similarity theory require a surface roughness to function. The system wide default is 1.0m, unless otherwise specified.

```@eval
using ..GasDispersion
using Plots

u0, z0 = 10, 10
zs = 0.1:.1:2
Businger = BasicEquationSet(BusingerWind(),nothing,nothing,nothing)

for (cls_lbl, class) in [("A",ClassA()),("B",ClassB()),("C",ClassC()),("D",ClassD()),("E",ClassE()),("F",ClassF())]
    p = plot(title="Windspeed for $class stability",ylabel="Relative Elevation (z/z_R)", xlabel="Relative Windspeed (u/u_R)")
    atm = SimpleAtmosphere(;windspeed=u0,windspeed_height=z0,stability=class)

    for row in [("Default",DefaultSet),("Irwin Rural",CCPSRural),("Irwin Urban",CCPSUrban),
                ("ISC3 Urban",ISC3Urban),("TNO",TNOPlume),("Businger et al.",Businger)]
        lbl, eqn = row
        u = [ GasDispersion.windspeed(atm,z*z0,eqn) for z in zs ]./u0
        plot!(p,u,zs,label=lbl)
    end

    savefig(p,"Class_$(cls_lbl)_windspeed.svg")

end

nothing
```

![](Class_A_windspeed.svg)
![](Class_B_windspeed.svg)
![](Class_C_windspeed.svg)
![](Class_D_windspeed.svg)
![](Class_E_windspeed.svg)
![](Class_F_windspeed.svg)

## Plume Dispersion

Plume dispersion parameters, $\sigma_y$ and $\sigma_z$ are functions of downwind distance and can take many different forms from simple power-law relations to complex piece-wise functions. The plume equation sets implement the plume dispersion parameters along with the windspeed correlations given above.

### Crosswind Dispersion

There are six correlations for the crosswind dispersion.
+ `Defaultσy` -- see [the defition in plume models](plume.md)
+ `BriggsRuralσy` -- [Briggs 1973](references.md), Appendix D
+ `BriggsUrbanσy` -- [Briggs 1973](references.md), Appendix D
+ `ISC3Ruralσy` -- [US EPA 1995](references.md), equation 1-32 and Table 1-1
+ `TNOPlumeσy` -- *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md)), Table 4-8
+ `Turnerσy` -- A set of digitized curves based on [Turner (1970)](references.md), as presented in [Lees 1996](references.md)

!!! note "Briggs Correlations"
    The Briggs correlations are in terms of $R$ and have been converted to $\sigma$ per [Griffiths (1994)](references.md).

```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for (cls_lbl, class) in [("A",ClassA()),("B",ClassB()),("C",ClassC()),("D",ClassD()),("E",ClassE()),("F",ClassF())]
    p = plot(title="Crosswind Dispersion for $class stability",ylabel="Crosswind Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for row in [("Default",DefaultSet),("Briggs Rural",CCPSRural),("Briggs Urban",CCPSUrban),
                ("ISC3 Rural",ISC3Rural),("TNO",TNOPlume),
                ("Turner",Turner)]
        lbl, eqn = row
        s = [ GasDispersion.crosswind_dispersion(x, class, eqn) for x in xs ]
        plot!(p,xs,s,label=lbl)
    end

    savefig(p,"Class_$(cls_lbl)_plume_sigma_y.svg")

end

nothing
```

![](Class_A_plume_sigma_y.svg)
![](Class_B_plume_sigma_y.svg)
![](Class_C_plume_sigma_y.svg)
![](Class_D_plume_sigma_y.svg)
![](Class_E_plume_sigma_y.svg)
![](Class_F_plume_sigma_y.svg)

### Vertical Dispersion

There are six correlations for the crosswind dispersion.
+ `Defaultσz` -- see [the defition in plume models](plume.md)
+ `BriggsRuralσz` -- [Briggs 1973](references.md), Appendix D
+ `BriggsUrbanσz` -- [Briggs 1973](references.md), Appendix D
+ `ISC3Ruralσz` -- [US EPA 1995](references.md), equation 1-34 and Table 1-2
+ `TNOPlumeσz` -- *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md)), Table 4-8 with a default surface roughnes of 0.1m
+ `Turnerσz` -- A set of digitized curves based on [Turner (1970)](references.md), as presented in [Lees 1996](references.md)

!!! note "Briggs Correlations"
    The Briggs correlations are in terms of $R$ and have been converted to $\sigma$ per [Griffiths (1994)](references.md).


```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for (cls_lbl, class) in [("A",ClassA()),("B",ClassB()),("C",ClassC()),("D",ClassD()),("E",ClassE()),("F",ClassF())]
    p = plot(title="Vertical Dispersion for $class stability",ylabel="Vertical Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for row in [("Default",DefaultSet),("Briggs Rural",CCPSRural),("Briggs Urban",CCPSUrban),
                ("ISC3 Rural",ISC3Rural),("TNO",TNOPlume),
                ("Turner",Turner)]
        lbl, eqn = row
        s = [ GasDispersion.vertical_dispersion(x, class, eqn) for x in xs ]
        plot!(p,xs,s,label=lbl)
    end

    savefig(p,"Class_$(cls_lbl)_plume_sigma_z.svg")

end

nothing
```

![](Class_A_plume_sigma_z.svg)
![](Class_B_plume_sigma_z.svg)
![](Class_C_plume_sigma_z.svg)
![](Class_D_plume_sigma_z.svg)
![](Class_E_plume_sigma_z.svg)
![](Class_F_plume_sigma_z.svg)


## Puff Dispersion
Puff dispersion parameters, $\sigma_x$, $\sigma_y$ and $\sigma_z$ are functions of the downwind distance to the cloud (puff) center and are generally given as power law relations. There are many fewer sources for these. The Puff equation sets implement these dispersion parameters along with the windspeed correlations given above.

### Downwind Dispersion

There two correlations for the downwind dispersion.
+ `CCPSPuffσx` -- [AIChE/CCPS 1999](references.md), Table 2-13
+ `TNOPuffσx` -- *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md)), pg 4.75


```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for (cls_lbl, class) in [("A",ClassA()),("B",ClassB()),("C",ClassC()),("D",ClassD()),("E",ClassE()),("F",ClassF())]
    p = plot(title="Downwind Dispersion for $class stability",ylabel="Downwind Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for row in [("CCPS",CCPSPuffRural),("TNO",TNOPuff)]
        lbl, eqn = row
        s = [ GasDispersion.downwind_dispersion(x, class, eqn) for x in xs ]
        plot!(p,xs,s,label=lbl)
    end

    savefig(p,"Class_$(cls_lbl)_puff_sigma_x.svg")

end

nothing
```

![](Class_A_puff_sigma_x.svg)
![](Class_B_puff_sigma_x.svg)
![](Class_C_puff_sigma_x.svg)
![](Class_D_puff_sigma_x.svg)
![](Class_E_puff_sigma_x.svg)
![](Class_F_puff_sigma_x.svg)

### Crosswind Dispersion

There two correlations for the crosswind dispersion.
+ `CCPSPuffσy` -- [AIChE/CCPS 1999](references.md), Table 2-13
+ `TNOPuffσy` -- *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md)), pg 4.75

```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for (cls_lbl, class) in [("A",ClassA()),("B",ClassB()),("C",ClassC()),("D",ClassD()),("E",ClassE()),("F",ClassF())]
    p = plot(title="Crosswind Dispersion for $class stability",ylabel="Crosswind Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for row in [("CCPS",CCPSPuffRural),("TNO",TNOPuff)]
        lbl, eqn = row
        s = [ GasDispersion.crosswind_dispersion(x, class, eqn) for x in xs ]
        plot!(p,xs,s,label=lbl)
    end

    savefig(p,"Class_$(cls_lbl)_puff_sigma_y.svg")

end

nothing
```

![](Class_A_puff_sigma_y.svg)
![](Class_B_puff_sigma_y.svg)
![](Class_C_puff_sigma_y.svg)
![](Class_D_puff_sigma_y.svg)
![](Class_E_puff_sigma_y.svg)
![](Class_F_puff_sigma_y.svg)

### Vertical Dispersion

There two correlations for the vertical dispersion.
+ `CCPSPuffσz` -- [AIChE/CCPS 1999](references.md), Table 2-13
+ `TNOPuffσz` -- *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md)), pg 4.75

Though in practice there are only two: the CCPS correlations do not distinguish between urban and rural locations for puff dispersion, and the default correlations *are* the CCPS correlations.

```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for (cls_lbl, class) in [("A",ClassA()),("B",ClassB()),("C",ClassC()),("D",ClassD()),("E",ClassE()),("F",ClassF())]
    p = plot(title="Vertical Dispersion for $class stability",ylabel="Vertical Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for row in [("CCPS",CCPSPuffRural),("TNO",TNOPuff)]
        lbl, eqn = row
        s = [ GasDispersion.vertical_dispersion(x, class, eqn) for x in xs ]
        plot!(p,xs,s,label=lbl)
    end

    savefig(p,"Class_$(cls_lbl)_puff_sigma_z.svg")

end

nothing
```

![](Class_A_puff_sigma_z.svg)
![](Class_B_puff_sigma_z.svg)
![](Class_C_puff_sigma_z.svg)
![](Class_D_puff_sigma_z.svg)
![](Class_E_puff_sigma_z.svg)
![](Class_F_puff_sigma_z.svg)