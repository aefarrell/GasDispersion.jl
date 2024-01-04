# Building Scenarios

A `Scenario` is a container for all of the information that a model may need to
produce a solution. The intention is for the `Scenario` to be re-usable, so that
the user may run the same scenario with multiple models without much difficulty.
Models also have specific parameters, those are handled in the model itself.

A `scenario_builder` function exists to help create valid `Scenario`s for
various standard release scenarios.
```@docs
scenario_builder

scenario_builder(::Substance, ::Type{JetSource}, ::Atmosphere)
```
