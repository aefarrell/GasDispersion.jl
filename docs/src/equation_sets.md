# Equation Sets

The various dispersion models each depend upon several parameters which are themselves, often, correlations. For any given parameter there are several different correlations in the literature. To make this more transparent, sets of correlations from standard texts have been prepared (in addition to the default correlations), allowing the user to *specify* which set to use.

Not all equation sets cover all of the possible correlations, unless otherwise specified the default correlations are used. Some models use internal models for the wind profile and dispersion, in such cases changing the equation set has no effect.

## Windspeed

For `SimpleAtmosphere`s the windspeed is determined using a power-law relationship:
```math
 u = u_{R} \left( z \over z_{R} \right)^{p}
```

There are five equation sets that give *p* as a function of stability class
+ `DefaultSet` -- see [the defition in release scenarios](scenarios.md)
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

```@eval
using ..GasDispersion
using Plots

u0, z0 = 1, 1
zs = 0:.1:10

for class in [ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]
    p = plot(title="Windspeed for $class stability",xlabel="Relative Elevation (z/z_R)", ylabel="Relative Windspeed (u/u_R)")

    for eqn in [DefaultSet,CCPSRural,CCPSUrban,ISC3Rural,ISC3Urban]
        u = [ GasDispersion._windspeed(u0,z0,z,class,eqn) for z in zs ]
        plot!(p,u,zs,label=eqn)
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

Plume dispersion parameters, $\sigma_y$ and $\sigma_z$ are functions of downwind distance and can take many different forms from simple power-law relations to complex piece-wise functions.

### Crosswind Dispersion

There are seven equation sets that set correlations for the crosswind dispersion.
+ `DefaultSet` -- see [the defition in plume models](plume.md)
+ `CCPSRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSUrban` -- [AIChE/CCPS 1999](references.md)
+ `ISC3Rural` -- [US EPA 1995](references.md)
+ `ISC3Urban` -- [US EPA 1995](references.md)
+ `TNO` -- These are the correlations given in the *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))
+ `Turner` -- A set of digitized curves based on [Turner (1970)](references.md), as presented in [Lees 1996](references.md)

```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for class in [ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]
    p = plot(title="Crosswind Dispersion for $class stability",ylabel="Crosswind Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for eqn in [DefaultSet,CCPSRural,CCPSUrban,ISC3Rural,ISC3Urban,TNO,Turner]
        s = [ GasDispersion.crosswind_dispersion(x, Plume, class, eqn) for x in xs ]
        plot!(p,xs,s,label=eqn)
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
+ `TNO` -- These are the correlations given in the *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))
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

    for eqn in [DefaultSet,CCPSRural,CCPSUrban,ISC3Rural,ISC3Urban,TNO,Turner]
        s = [ GasDispersion.vertical_dispersion(x, Plume, class, eqn) for x in xs ]
        plot!(p,xs,s,label=eqn)
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
Puff dispersion parameters, $\sigma_x$, $\sigma_y$ and $\sigma_z$ are functions of the downwind distance to the cloud (puff) center and are generally given as power law relations. There are many fewer sources for these.

### Downwind Dispersion

There are four equation sets that set correlations for the downwind dispersion.
+ `DefaultSet` -- see [the defition in plume models](plume.md)
+ `CCPSRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSUrban` -- [AIChE/CCPS 1999](references.md)
+ `TNO` -- These are the correlations given in the *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))

Though in practice there are only two: the CCPS correlations do not distinguish between urban and rural locations for puff dispersion, and the default correlations *are* the CCPS correlations.

```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for class in [ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]
    p = plot(title="Downwind Dispersion for $class stability",ylabel="Downwind Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for eqn in [DefaultSet,CCPSRural,CCPSUrban,TNO]
        s = [ GasDispersion.downwind_dispersion(x, Puff, class, eqn) for x in xs ]
        plot!(p,xs,s,label=eqn)
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
+ `DefaultSet` -- see [the defition in plume models](plume.md)
+ `CCPSRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSUrban` -- [AIChE/CCPS 1999](references.md)
+ `TNO` -- These are the correlations given in the *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))

Though in practice there are only two: the CCPS correlations do not distinguish between urban and rural locations for puff dispersion, and the default correlations *are* the CCPS correlations.

```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for class in [ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]
    p = plot(title="Crosswind Dispersion for $class stability",ylabel="Crosswind Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for eqn in [DefaultSet,CCPSRural,CCPSUrban,TNO]
        s = [ GasDispersion.crosswind_dispersion(x, Puff, class, eqn) for x in xs ]
        plot!(p,xs,s,label=eqn)
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
+ `DefaultSet` -- see [the defition in plume models](plume.md)
+ `CCPSRural` -- [AIChE/CCPS 1999](references.md)
+ `CCPSUrban` -- [AIChE/CCPS 1999](references.md)
+ `TNO` -- These are the correlations given in the *TNO Yellow Book* ([Bakkum and Duijm 2005](references.md))

Though in practice there are only two: the CCPS correlations do not distinguish between urban and rural locations for puff dispersion, and the default correlations *are* the CCPS correlations.

```@eval
using ..GasDispersion
using Plots

xs = range(1,1e5,1000)

for class in [ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]
    p = plot(title="Vertical Dispersion for $class stability",ylabel="Vertical Dispersion, m", xlabel="Downwind distance, m", xaxis=:log10, yaxis=:log10, legend=:bottomright)

    for eqn in [DefaultSet,CCPSRural,CCPSUrban,TNO]
        s = [ GasDispersion.vertical_dispersion(x, Puff, class, eqn) for x in xs ]
        plot!(p,xs,s,label=eqn)
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