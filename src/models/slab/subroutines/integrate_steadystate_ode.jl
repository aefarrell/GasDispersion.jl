"""OrdinaryDiffEq backend for the steady-state phase."""

"""
Initialize the full-domain solve context.
"""
function _slab_steady_integrator(solver::OrdinaryDiffEqAlgorithm,
                                 input::SLAB_Steady_IntegratorInput)
    controls = input.controls
    context = _slab_steady_context(input.params, input.base, input.idpf, nothing;
        rmi=controls.rmi, alfg=controls.alfg, sru0=controls.sru0, bbx=controls.bbx,
        x0=input.x0)
    return OrdinaryDiffEqIntegratorContext(solver, context)
end

"""
Integrate across the full spatial domain with adaptive timesteps.

The callback rebases the reduced ODE state after each accepted step so the
thermodynamic state and controls used by the RHS stay current. Output samples
are interpolated on a logarithmic spatial grid, which does not constrain the
solver's accepted timesteps.
"""
function _slab_steady_integrate!(intctx::OrdinaryDiffEqIntegratorContext,
        vecs, vars, params, input::SLAB_Steady_IntegratorInput;
        nxtr, cv, ug, tim, bbx, bx, betax)
    context = intctx.context
    F = eltype(vecs.x)
    targets = Vector{F}(undef, vars.mffm - vars.nxi + 1)
    xf = input.x1
    targets .= exp.(range(log(input.x0), log(xf); length=length(targets)+1))[2:end]

    segments = SLAB_Steady_ODESegment{F,typeof(context)}[]
    last_t = Ref(input.x0)
    last_dt = Ref(zero(F))
    terminal_state = Ref{Union{Nothing,SLAB_Steady_Phase_State{F}}}(nothing)
    terminal_cv = Ref(zero(F))
    reference = Ref(_slab_steady_reference_state(input.base))
    controls = Ref(input.controls)
    affect! = function (integrator)
        t = integrator.t
        if t > last_t[]
            push!(segments, SLAB_Steady_ODESegment(last_t[], t,
                _slab_steady_context_snapshot(context)))
            last_dt[] = t - last_t[]
        end
        state, derived = _slab_steady_project(integrator.u, context, t)
        reference[], controls[] = _slab_steady_reference_update(
            SLAB_Steady_StepResult(state, derived.cv), controls[], params.met.rhoa)
        context.base = _slab_steady_loop_state(state)
        context.y0 = integrator.u
        context.x0 = t
        context.rmi = controls[].rmi
        context.alfg = controls[].alfg
        context.sru0 = controls[].sru0
        context.bbx = controls[].bbx
        last_t[] = t
        if state.qint >= 0.5*params.spl.qtcs
            terminal_state[] = context.base
            terminal_cv[] = derived.cv
            terminate!(integrator)
        end
        return nothing
    end
    callback = DiscreteCallback((u,t,integrator) -> true, affect!;
        save_positions=(false,false))
    problem = ODEProblem(_slab_steady_rhs, context.y0, (input.x0, xf), context)
    solution = solve(problem, intctx.solver;
        merge(input.solver_kwargs, (;callback=callback, dense=true))...)

    stop_t = solution.t[end]
    if last_t[] < stop_t
        push!(segments, SLAB_Steady_ODESegment(last_t[], stop_t,
            _slab_steady_context_snapshot(context)))
    end

    _state_at(x) = begin
        index = findlast(segment -> segment.x0 <= x <= segment.x1, segments)
        index === nothing && error("No ODE context was saved for steady-state position $x")
        state, derived = _slab_steady_project(solution(x), segments[index].context, x)
        return state, derived.cv
    end
    stored_nxtr = nxtr
    for (offset, sample_x) in enumerate(targets)
        nx = vars.nxi + offset - 1
        if terminal_state[] !== nothing && sample_x >= stop_t && stored_nxtr == nxtr
            state = terminal_state[]
            _slab_sub_store!(vecs,nx,stop_t,state.bb,state.b,state.vg,state.cm,state.t,
                state.rho,state.u,state.h,terminal_cv[],state.beta,state.w,state.v,
                params.met.cmdaa,state.cmw,state.cmwv,state.cmev,state.uab,state.wc,
                state.zc,state.qint,tim,bbx,bx,betax,ug,state.vx)
            vecs.tccp[nx] = (state.qint+state.qint)/params.spl.qs
            stored_nxtr = nx
            break
        elseif sample_x <= stop_t
            state, sample_cv = _state_at(sample_x)
            _slab_sub_store!(vecs,nx,sample_x,state.bb,state.b,state.vg,state.cm,
                state.t,state.rho,state.u,state.h,sample_cv,state.beta,state.w,state.v,
                params.met.cmdaa,state.cmw,state.cmwv,state.cmev,state.uab,state.wc,
                state.zc,state.qint,tim,bbx,bx,betax,ug,state.vx)
            vecs.tccp[nx] = (state.qint+state.qint)/params.spl.qs
        else
            break
        end
    end

    final_state = terminal_state[] === nothing ? context.base : terminal_state[]
    final_cv = terminal_state[] === nothing ? cv : terminal_cv[]
    return (base=final_state, x=stop_t, dx=last_dt[], cv=final_cv,
        nxtr=stored_nxtr, reference=reference[], controls=controls[],
        stopped=terminal_state[] !== nothing, ode_solution=solution,
        ode_segments=segments)
end

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