"""Persistent OrdinaryDiffEq driver for the steady-state phase."""

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


function _slab_steady_ode_segment(intctx::OrdinaryDiffEqIntegratorContext, x0, x1)
    return SLAB_Steady_ODESegment(x0, x1, deepcopy(intctx.context))
end

function _slab_steady_ode_segments(intctx::OrdinaryDiffEqIntegratorContext,
                                   ::Type{F}) where {F}
    return SLAB_Steady_ODESegment{F,typeof(intctx.context)}[]
end

_slab_steady_solution(intctx::OrdinaryDiffEqIntegratorContext) = intctx.integrator.sol