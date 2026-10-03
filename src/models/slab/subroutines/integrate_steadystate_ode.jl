"""
OrdinaryDiffEq backend for the steady-state phase.

One persistent ODE integrator is reused across the shared loop's substeps.
Each substep reinitializes its independent variables and updates the mutable
RHS context; saved segment contexts preserve the coefficients needed for later
field interpolation.
"""

"""
Initialize an OrdinaryDiffEq integrator for the first steady substep.

The initial context and ODE problem are reused for all later substeps; the
integrator is configured here with the caller-provided solver options.
"""
function _slab_steady_integrator(solver::OrdinaryDiffEqAlgorithm,
                                 input::SLAB_Steady_IntegratorInput)
    controls = input.controls
    context = _slab_steady_context(input.params, input.base, input.idpf, nothing;
        rmi=controls.rmi, alfg=controls.alfg, sru0=controls.sru0, bbx=controls.bbx,
        x0=input.x0)
    problem = ODEProblem(_slab_steady_rhs, context.y0, (input.x0, input.x1), context)
    integrator = init(problem, solver; input.solver_kwargs...)
    return OrdinaryDiffEqIntegratorContext(integrator,context)
end

"""
Advance the persistent ODE integrator over one requested spatial interval.

Update the RHS context from the current step input, reinitialize the solver at
the interval's starting point, integrate through its endpoint, then project the
final reduced ODE state into the full SLAB phase state.
"""
function _slab_steady_step!(intctx::OrdinaryDiffEqIntegratorContext,
                            step::SLAB_Steady_StepInput)
    input = step.input
    base = input.base
    controls = input.controls
    x, xf = input.x0, input.x1
    integrator, context = intctx.integrator, intctx.context
    context.base = base
    context.x0 = x
    context.rmi = controls.rmi
    context.alfg = controls.alfg
    context.sru0 = controls.sru0
    context.bbx = controls.bbx
    context.y0 = SVector{11,typeof(base.r)}(
        base.r, base.bbv, base.bv, base.g, base.sft, base.sfu,
        base.gw, base.zc, base.qint, base.sfy, base.sfz)
    reinit!(integrator, context.y0; t0=x, tf=xf, erase_sol=false)
    requested_dt = get(input.solver_kwargs, :dt, xf - x)
    set_proposed_dt!(integrator, min(requested_dt, xf - x))
    while integrator.t < xf
        step!(integrator)
    end
    state, derived = _slab_steady_project(integrator.u, context, xf)
    return SLAB_Steady_StepResult(state, derived.cv)
end

"""
Save the RHS context associated with one completed ODE interval.

The integrator mutates its live context on the next step, so each segment needs
a snapshot. The snapshot copies the small context object but shares immutable
parameters and their arrays instead of deep-copying them for every interval.
"""
function _slab_steady_append_ode_segment!(segments,
        intctx::OrdinaryDiffEqIntegratorContext, x0, x1)
    push!(segments, SLAB_Steady_ODESegment(
        x0, x1, _slab_steady_context_snapshot(intctx.context)))
    return nothing
end

"""Allocate the segment collection with the context type used by this integrator."""
function _slab_steady_ode_segments(intctx::OrdinaryDiffEqIntegratorContext,
                                   ::Type{F}) where {F}
    return SLAB_Steady_ODESegment{F,typeof(intctx.context)}[]
end

_slab_steady_solution(intctx::OrdinaryDiffEqIntegratorContext) = intctx.integrator.sol

"""
Evaluate an ODE-backed cloud-field interpolation at distance `x`.

Use the saved segment context to reconstruct the full thermodynamic state from
the dense ODE solution, then derive the requested cloud field. Outside the
integrated segments, delegate to the tabulated interpolation fallback.
"""
function (itp::SLAB_ODE_FieldInterpolation)(x)
    index = findlast(segment -> segment.x0 <= x <= segment.x1, itp.segments)
    index === nothing && return itp.fallback(x)
    y = itp.solution(x)
    state, _ = _slab_steady_project(y, itp.segments[index].context, x)
    if itp.field === :b
        return state.b
    elseif itp.field === :zc
        return state.zc
    end
    p = itp.params
    cv = (p.met.wmae*state.cm)/(p.rgp.wms+(p.met.wmae-p.rgp.wms)*state.cm)
    point = _slab_editcc_point(p,x,itp.x0,state.zc,state.h,state.b,state.beta,
        state.uab,state.cm,cv,itp.bx_x(x),itp.bbx_x(x),itp.tcld(x))
    return getproperty(point,itp.field)
end

"""
Construct output interpolations for an OrdinaryDiffEq steady solution.

The stored SLAB vectors provide the compatibility fallback and fields not
reconstructed from the ODE. For integrated intervals, dense-solution wrappers
provide the cloud fields evaluated from the segment-specific RHS contexts.
"""
function _slab_output_interpolations(cc::SLAB_CC_Vecs,
        steady_solution::SLAB_Steady_Solution{I,F,A,O,G}, params) where {
        I,F,A,O,G<:AbstractVector{<:SLAB_Steady_ODESegment}}
    legacy = _slab_legacy_interpolations(cc)
    if isempty(steady_solution.ode_segments)
        return legacy
    end
    makefield(field, fallback) = SLAB_ODE_FieldInterpolation(
        steady_solution.ode_solution, steady_solution.ode_segments, params, field,
        cc.x[1], legacy.bx_x, legacy.bbx_x, legacy.tcld, fallback)
    return SLAB_Interpolations(
        makefield(:cc,legacy.cc), makefield(:b,legacy.b),
        makefield(:betac,legacy.betac), makefield(:zc,legacy.zc),
        makefield(:sig,legacy.sig), legacy.xc, legacy.bx, legacy.betax,
        legacy.bx_x, legacy.bbx_x, legacy.tcld, steady_solution.ode_solution)
end