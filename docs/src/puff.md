# Puff Models

Puff models are for "instantaneous" releases or other time-dependent releases.

```@docs
puff
```

## Gaussian Puffs

A Gaussian puff model is defined by gaussian distributions in the x, y, and z
directions and travels downwind at the ambient windspeed. The dispersion
parameters of the gaussians are correlated with the downwind distance of cloud
center.

```@docs
puff(::Scenario, ::Type{GaussianPuff})
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