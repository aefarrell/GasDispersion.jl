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

There are four equation sets for puff models which define the correlations for windspeed, downwind dispersion, crosswind dispersion, and vertical dispersion:
+ `DefaultPuffSet` -- see [the defition in plume models](plume.md)
+ `CCPSPuffRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSPuffUrban` -- [AIChE/CCPS 1999](references.md)
+ `TNOPuff` -- *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))

An equation set intended for puff models can be used for a plume model, but not vice-versa unless otherwise noted. Equation sets for plume models do not typically define a downwind dispersion, $\sigma_x$, and without some additional details on how to handle that there is no way to calculate puff dispersion.

It is possible to define one's own equation set by either mixing and matching existing correlations or defining singletons for the new correlations and overloading the appropriate methods internal to `GasDispersion`, but this is dangerous as the internals are subject to change. If it is a common equation set please feel free to add it to `GasDispersion` by initiating a pull request.

### Example Usage

Using the same example scenario for the basic gaussian plume, suppose we want to see how this changes as we change correlations. Going back to the scenario definition:

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

The plume using the default equation set is simply this

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
dflt = plume(scn, GaussianPlume, DefaultSet())

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.NoPlumeRise, ClassF, GasDispersion.BasicEquationSet{GasDispersion.DefaultWind, Nothing, GasDispersion.Defaultσy, GasDispersion.Defaultσz}}(Scenario{Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF)), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 1.150112899011524, 3.5, GasDispersion.NoPlumeRise(), ClassF, GasDispersion.BasicEquationSet{GasDispersion.DefaultWind, Nothing, GasDispersion.Defaultσy, GasDispersion.Defaultσz}())

```

For each of the plume equation sets we can easily create corresponding plume solutions

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
ccps_rurl = plume(scn, GaussianPlume, CCPSRural())

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.NoPlumeRise, ClassF, GasDispersion.BasicEquationSet{GasDispersion.IrwinRural, Nothing, GasDispersion.BriggsRuralσy, GasDispersion.BriggsRuralσz}}(Scenario{Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF)), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 0.8420321686971456, 3.5, GasDispersion.NoPlumeRise(), ClassF, GasDispersion.BasicEquationSet{GasDispersion.IrwinRural, Nothing, GasDispersion.BriggsRuralσy, GasDispersion.BriggsRuralσz}())

```

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
ccps_urb = plume(scn, GaussianPlume, CCPSUrban())

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.NoPlumeRise, ClassF, GasDispersion.BasicEquationSet{GasDispersion.IrwinUrban, Nothing, GasDispersion.BriggsUrbanσy, GasDispersion.BriggsUrbanσz}}(Scenario{Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF)), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 0.7989729675905327, 3.5, GasDispersion.NoPlumeRise(), ClassF, GasDispersion.BasicEquationSet{GasDispersion.IrwinUrban, Nothing, GasDispersion.BriggsUrbanσy, GasDispersion.BriggsUrbanσz}())

```

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
isc3_rurl = plume(scn, GaussianPlume, ISC3Rural())

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.NoPlumeRise, ClassF, GasDispersion.BasicEquationSet{GasDispersion.IrwinRural, Nothing, GasDispersion.ISC3Ruralσy, GasDispersion.ISC3Ruralσz}}(Scenario{Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF)), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 0.8420321686971456, 3.5, GasDispersion.NoPlumeRise(), ClassF, GasDispersion.BasicEquationSet{GasDispersion.IrwinRural, Nothing, GasDispersion.ISC3Ruralσy, GasDispersion.ISC3Ruralσz}())

```

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
isc3_urb = plume(scn, GaussianPlume, ISC3Urban())

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.NoPlumeRise, ClassF, GasDispersion.BasicEquationSet{GasDispersion.ISC3UrbanWind, Nothing, GasDispersion.BriggsUrbanσy, GasDispersion.BriggsUrbanσz}}(Scenario{Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF)), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 1.0947417281650496, 3.5, GasDispersion.NoPlumeRise(), ClassF, GasDispersion.BasicEquationSet{GasDispersion.ISC3UrbanWind, Nothing, GasDispersion.BriggsUrbanσy, GasDispersion.BriggsUrbanσz}())

```

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
tno = plume(scn, GaussianPlume, TNOPlume())

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.NoPlumeRise, ClassF, GasDispersion.BasicEquationSet{GasDispersion.TNOWind, Nothing, GasDispersion.TNOPlumeσy, GasDispersion.TNOPlumeσz}}(Scenario{Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF)), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 0.6446987024117233, 3.5, GasDispersion.NoPlumeRise(), ClassF, GasDispersion.BasicEquationSet{GasDispersion.TNOWind, Nothing, GasDispersion.TNOPlumeσy, GasDispersion.TNOPlumeσz}())

```

```jldoctest eqnset_example; output = false, filter = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
turner = plume(scn, GaussianPlume, Turner())

# output

GasDispersion.GaussianPlumeSolution{Float64, GasDispersion.NoPlumeRise, ClassF, GasDispersion.BasicEquationSet{GasDispersion.DefaultWind, Nothing, GasDispersion.Turnerσy, GasDispersion.Turnerσz}}(Scenario{Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}, HorizontalJet{Float64}, SimpleAtmosphere{Float64, ClassF}}(Substance{String, GasDispersion.Antoine{Float64}, Float64, Float64, Float64, Int64, Int64, Int64}("propane", 0.044096, GasDispersion.Antoine{Float64}(9.773719865868816, 2257.9247634130143, 0.0), 1.864931992847327, 526.13, 288.15, 101325.0, 1.142, 231.02, 425740, 1678, 2520), HorizontalJet{Float64}(0.08991798763471508, Inf, 0.01, 208.10961399327573, 3.5, 288765.2212333958, 278.3846872082166, 0.0), SimpleAtmosphere{Float64, ClassF}(101325.0, 298.15, 1.5, 10.0, 0.0, ClassF)), :gaussian, 0.08991798763471508, 0.9999999999999998, 1.8023818673116125, 1.150112899011524, 3.5, GasDispersion.NoPlumeRise(), ClassF, GasDispersion.BasicEquationSet{GasDispersion.DefaultWind, Nothing, GasDispersion.Turnerσy, GasDispersion.Turnerσz}())

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

scn = scenario_builder(propane, JetSource; 
       phase = :gas,
       diameter = 0.01,  # m
       dischargecoef = 0.85,
       temperature = T1, # K
       pressure = P1,    # Pa
       height = 3.5)     # m, height of hole above the ground

dflt = plume(scn, GaussianPlume, DefaultSet())
ccps_rurl = plume(scn, GaussianPlume, CCPSRural())
ccps_urb = plume(scn, GaussianPlume, CCPSUrban())
isc3_rurl = plume(scn, GaussianPlume, ISC3Rural())
isc3_urb = plume(scn, GaussianPlume, ISC3Urban())
turner = plume(scn, GaussianPlume, Turner())

solset = [("Default", dflt), ("CCPSRural", ccps_rurl), 
          ("CCPSUrban", ccps_urb), ("ISC3Rural", isc3_rurl),
          ("ISC3Urban", isc3_urb), ("Turner", turner)]

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

There are five equation sets that give *p* as a function of stability class
+ `DefaultSet` -- see [the defition in release scenarios](scenarios.md), also what is used for `Turner`
+ `CCPSRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSUrban` -- [AIChE/CCPS 1999](references.md)
+ `ISC3Rural` -- [US EPA 1995](references.md)
+ `ISC3Urban` -- [US EPA 1995](references.md)

| Stability Class | `DefaultSet`  | `CCPSRural` | `CCPSUrban` | `ISC3Rural` | `ISC3Urban` |
|:---------------:|:-------------:|:-----------:|:-----------:|:-----------:|:-----------:|
|        A        |  0.108        |  0.07       |  0.15       |  0.07       |  0.15       |
|        B        |  0.112        |  0.07       |  0.15       |  0.07       |  0.15       |
|        C        |  0.120        |  0.10       |  0.20       |  0.10       |  0.20       |
|        D        |  0.142        |  0.15       |  0.25       |  0.15       |  0.25       |
|        E        |  0.203        |  0.35       |  0.40       |  0.35       |  0.30       |
|        F        |  0.253        |  0.55       |  0.60       |  0.55       |  0.30       |

!!! note "Note"
    The CCPS and ISC3 correlations both use the same parameters for Rural settings but differ for Urban settings with stable atmospheres (class E and F) 

There is only one equation set which uses a logarithmic profile based on Monin-Obukhov similarity theory.
+ `TNOPlume` -- *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))

!!! note "Note"
    The logarithmic profiles, such as used in `TNOPlume`, can produce bizarre answers for plumes close to the ground. They also struggle with windspeeds measured at elevations very near the ground.

```@eval
using ..GasDispersion
using Plots

u0, z0 = 10, 10
zs = 0.1:.1:2

for class in [ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]
    p = plot(title="Windspeed for $class stability",ylabel="Relative Elevation (z/z_R)", xlabel="Relative Windspeed (u/u_R)")
    atm = SimpleAtmosphere(;windspeed=u0,windspeed_height=z0,stability=class)

    for row in [("DefaultSet",DefaultSet()),("CCPSRural",CCPSRural()),("CCPSUrban",CCPSUrban()),
                ("ISC3Rural",ISC3Rural()),("ISC3Urban",ISC3Urban()),("TNOPlume",TNOPlume())]
        lbl, eqn = row
        u = [ GasDispersion._windspeed(atm,z*z0,eqn) for z in zs ]./u0
        plot!(p,u,zs,label=lbl)
    end

    savefig(p,"$(class)_windspeed.svg")

end

nothing
```

![](ClassA_windspeed.svg)
![](ClassB_windspeed.svg)
![](ClassC_windspeed.svg)
![](ClassD_windspeed.svg)
![](ClassE_windspeed.svg)
![](ClassF_windspeed.svg)

## Plume Dispersion

Plume dispersion parameters, $\sigma_y$ and $\sigma_z$ are functions of downwind distance and can take many different forms from simple power-law relations to complex piece-wise functions. The plume equation sets implement the plume dispersion parameters along with the windspeed correlations given above.

### Crosswind Dispersion

There are seven equation sets that set correlations for the crosswind dispersion.
+ `DefaultSet` -- see [the defition in plume models](plume.md)
+ `CCPSRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSUrban` -- [AIChE/CCPS 1999](references.md)
+ `ISC3Rural` -- [US EPA 1995](references.md)
+ `ISC3Urban` -- [US EPA 1995](references.md)
+ `TNOPlume` -- These are the plume correlations given in the *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))
+ `Turner` -- A set of digitized curves based on [Turner (1970)](references.md), as presented in [Lees 1996](references.md)

```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for class in [ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]
    p = plot(title="Crosswind Dispersion for $class stability",ylabel="Crosswind Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for row in [("DefaultSet",DefaultSet()),("CCPSRural",CCPSRural()),("CCPSUrban",CCPSUrban()),
                ("ISC3Rural",ISC3Rural()),("ISC3Urban",ISC3Urban()),("TNOPlume",TNOPlume()),
                ("Turner",Turner())]
        lbl, eqn = row
        s = [ GasDispersion.crosswind_dispersion(x, class, eqn) for x in xs ]
        plot!(p,xs,s,label=lbl)
    end

    savefig(p,"$(class)_plume_sigma_y.svg")

end

nothing
```

![](ClassA_plume_sigma_y.svg)
![](ClassB_plume_sigma_y.svg)
![](ClassC_plume_sigma_y.svg)
![](ClassD_plume_sigma_y.svg)
![](ClassE_plume_sigma_y.svg)
![](ClassF_plume_sigma_y.svg)

### Vertical Dispersion

There are seven equation sets that set correlations for the crosswind dispersion.
+ `DefaultSet` -- see [the defition in plume models](plume.md)
+ `CCPSRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSUrban` -- [AIChE/CCPS 1999](references.md)
+ `ISC3Rural` -- [US EPA 1995](references.md)
+ `ISC3Urban` -- [US EPA 1995](references.md)
+ `TNOPlume` -- These are the plume correlations given in the *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))
+ `Turner` -- A set of digitized curves based on [Turner (1970)](references.md), as presented in [Lees 1996](references.md)

!!! note "Corrections"
    The CCPS correlations for the vertical plume dispersion, $\sigma_z$, in urban terrain 
    contains two typos. These have been corrected as per [Griffiths (1994)](references.md).


```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for class in [ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]
    p = plot(title="Vertical Dispersion for $class stability",ylabel="Vertical Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for row in [("DefaultSet",DefaultSet()),("CCPSRural",CCPSRural()),("CCPSUrban",CCPSUrban()),
                ("ISC3Rural",ISC3Rural()),("ISC3Urban",ISC3Urban()),("TNOPlume",TNOPlume()),
                ("Turner",Turner())]
        lbl, eqn = row
        s = [ GasDispersion.vertical_dispersion(x, class, eqn) for x in xs ]
        plot!(p,xs,s,label=lbl)
    end

    savefig(p,"$(class)_plume_sigma_z.svg")

end

nothing
```

![](ClassA_plume_sigma_z.svg)
![](ClassB_plume_sigma_z.svg)
![](ClassC_plume_sigma_z.svg)
![](ClassD_plume_sigma_z.svg)
![](ClassE_plume_sigma_z.svg)
![](ClassF_plume_sigma_z.svg)


## Puff Dispersion
Puff dispersion parameters, $\sigma_x$, $\sigma_y$ and $\sigma_z$ are functions of the downwind distance to the cloud (puff) center and are generally given as power law relations. There are many fewer sources for these. The Puff equation sets implement these dispersion parameters along with the windspeed correlations given above.

### Downwind Dispersion

There are four equation sets that set correlations for the downwind dispersion.
+ `DefaultPuffSet` -- see [the defition in plume models](plume.md)
+ `CCPSPuffRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSPuffUrban` -- [AIChE/CCPS 1999](references.md)
+ `TNOPuff` -- These are the correlations given in the *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))

Though in practice there are only two: the CCPS correlations do not distinguish between urban and rural locations for puff dispersion, and the default correlations *are* the CCPS correlations.

```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for class in [ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]
    p = plot(title="Downwind Dispersion for $class stability",ylabel="Downwind Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for row in [("DefaultPuffSet",DefaultPuffSet()),("CCPSPuffRural",CCPSPuffRural()),("CCPSPuffUrban",CCPSPuffUrban()),
                ("TNOPuff",TNOPuff())]
        lbl, eqn = row
        s = [ GasDispersion.downwind_dispersion(x, class, eqn) for x in xs ]
        plot!(p,xs,s,label=lbl)
    end

    savefig(p,"$(class)_puff_sigma_x.svg")

end

nothing
```

![](ClassA_puff_sigma_x.svg)
![](ClassB_puff_sigma_x.svg)
![](ClassC_puff_sigma_x.svg)
![](ClassD_puff_sigma_x.svg)
![](ClassE_puff_sigma_x.svg)
![](ClassF_puff_sigma_x.svg)

### Crosswind Dispersion

There are four equation sets that set correlations for the crosswind dispersion.
+ `DefaultPuffSet` -- see [the defition in plume models](plume.md)
+ `CCPSPuffRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSPuffUrban` -- [AIChE/CCPS 1999](references.md)
+ `TNOPuff` -- These are the correlations given in the *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))

Though in practice there are only two: the CCPS correlations do not distinguish between urban and rural locations for puff dispersion, and the default correlations *are* the CCPS correlations.

```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for class in [ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]
    p = plot(title="Crosswind Dispersion for $class stability",ylabel="Crosswind Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for row in [("DefaultPuffSet",DefaultPuffSet()),("CCPSPuffRural",CCPSPuffRural()),("CCPSPuffUrban",CCPSPuffUrban()),
                ("TNOPuff",TNOPuff())]
        lbl, eqn = row
        s = [ GasDispersion.crosswind_dispersion(x, class, eqn) for x in xs ]
        plot!(p,xs,s,label=lbl)
    end

    savefig(p,"$(class)_puff_sigma_y.svg")

end

nothing
```

![](ClassA_puff_sigma_y.svg)
![](ClassB_puff_sigma_y.svg)
![](ClassC_puff_sigma_y.svg)
![](ClassD_puff_sigma_y.svg)
![](ClassE_puff_sigma_y.svg)
![](ClassF_puff_sigma_y.svg)

### Vertical Dispersion

There are four equation sets that set correlations for the vertical dispersion.
+ `DefaultPuffSet` -- see [the defition in plume models](plume.md)
+ `CCPSPuffRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSPuffUrban` -- [AIChE/CCPS 1999](references.md)
+ `TNOPuff` -- These are the correlations given in the *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))

Though in practice there are only two: the CCPS correlations do not distinguish between urban and rural locations for puff dispersion, and the default correlations *are* the CCPS correlations.

```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for class in [ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]
    p = plot(title="Vertical Dispersion for $class stability",ylabel="Vertical Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for row in [("DefaultPuffSet",DefaultPuffSet()),("CCPSPuffRural",CCPSPuffRural()),("CCPSPuffUrban",CCPSPuffUrban()),
                ("TNOPuff",TNOPuff())]
        lbl, eqn = row
        s = [ GasDispersion.vertical_dispersion(x, class, eqn) for x in xs ]
        plot!(p,xs,s,label=lbl)
    end

    savefig(p,"$(class)_puff_sigma_z.svg")

end

nothing
```

![](ClassA_puff_sigma_z.svg)
![](ClassB_puff_sigma_z.svg)
![](ClassC_puff_sigma_z.svg)
![](ClassD_puff_sigma_z.svg)
![](ClassE_puff_sigma_z.svg)
![](ClassF_puff_sigma_z.svg)