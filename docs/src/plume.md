# Plume Models

Plume models are models of continuous, steady-state, releases and are time
independent, this includes, for example, emissions from elevated stacks.

```@docs
plume
```

## Gaussian Plumes

A Gaussian plume is a steady state plume defined by a Gaussian distribution in
the y and z directions and an exponential decay in the x direction. The
dispersion parameters are correlated to the downwind distance.

```@docs
plume(::Scenario, ::Type{GaussianPlume})
```

## Simple Jet Plumes

The simple jet plume is a steady state turbulent jet defined by a Gaussian
distribution in the y and z directions. It is similar to the Gaussian plume,
however in this case the momentum forming the plume comes entirely from the jet
and not the ambient windspeed.

```@docs
plume(::Scenario, ::Type{SimpleJet})
```

## Britter-McQuaid Model

The Britter-McQuaid model is an empirical correlation for dense plume
dispersion. The model generates an interpolation function for the centerline
concentration at the downwind distance x.


```@docs
plume(::Scenario, ::Type{BritterMcQuaidPlume})
```
