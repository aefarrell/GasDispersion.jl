# SLAB Follow-up

## ODE RK4 parity: INPR4 vertical jet

The fixed-step OrdinaryDiffEq `RK4` backend currently fails one assertion in the `INPR4 Vertical Jet` regression (`test/models/slab_tests.jl`). In the observed run, the largest mismatch was about `0.124 K` in temperature at output row 15; the other 14 of 15 SLAB assertions passed. The legacy solver remains the reference and default.

Trace the first divergent RK stage around the steady-to-transient handoff, comparing the legacy loop and ODE backend state, RHS, and handoff values. Preserve the reference test tolerance; add a focused regression for the identified state/transition behavior. Do not promote the ODE backend as equivalent until the INPR4 comparison passes.

## Refactor the transient phase

Apply the steady-phase architecture to the transient SLAB model:

- Represent the complete transient state explicitly and make the RHS a function of that state, parameters, and time.
- Keep the legacy transient RK implementation as the reference solver, dispatched separately from the OrdinaryDiffEq implementation.
- Use a persistent OrdinaryDiffEq integrator with manual stepping and expose solver algorithm/options through the existing SLAB solver interface.
- Retain the transient `ODESolution` and use its interpolation for ODE-backed transient output; retain Akima interpolation for legacy output.
- Validate the steady-to-transient initialization and compare all transient reference cases before enabling the ODE backend by default.

## Fix compat issues

Julia 1.3 no longer compiles the package. There is a compatibility issue with SciMLBase and OrdinaryDiffEqLowOrderRK. The real solution is to no longer support Julia 1.3 in future versions of GasDispersion, and find the minimal version of Julia which supports this new dependancy
