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

function _save_function(u, t, integrator)
    return _slab_steady_project(u, integrator.p, t)
end


"""
Integrate across the full spatial domain with adaptive timesteps.

The saving callback records each accepted state directly. A second callback
rebases the RHS context and terminates at the steady-to-transient threshold.
"""
function _slab_steady_integrate!(intctx::OrdinaryDiffEqIntegratorContext,
        vecs, vars, params, input::SLAB_Steady_IntegratorInput;
        ug, tim, bbx, bx, betax)
    context = intctx.context
    F = eltype(vecs.x)

    # The saving callback calculates the current state and saves it
    saved_values = SavedValues(F, SLAB_Steady_Phase_State{F})
    save_callback = SavingCallback(_save_function, saved_values;
        save_everystep=true, save_start=true, save_end=true)

    # A discrete call back refreshs the local context after each step
    # and checks whether the termination criteria has been met
    last_t = Ref(input.x0)
    last_dt = Ref(zero(F))
    terminal_state = Ref{Union{Nothing,SLAB_Steady_Phase_State{F}}}(nothing)
    reference = Ref(_slab_steady_reference_state(input.base))
    controls = Ref(input.controls)
    
    function update_state!(integrator)
        t = integrator.t
        if t > last_t[]
            last_dt[] = t - last_t[]
        end

        # Refresh derived SLAB state and the controls for the next RHS
        # evaluation, then stop once the transient threshold is met.
        state = _slab_steady_project(integrator.u, context, t)
        reference[], controls[] = _slab_steady_reference_update( state, controls[], params.met.rhoa)
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
            terminate!(integrator)
        end
        return nothing
    end

    update_callback = DiscreteCallback((u,t,integrator) -> true, update_state!;
        save_positions=(false,false))

    # Integrate the ODE problem 
    problem = ODEProblem(_slab_steady_rhs, context.y0,
                          (input.x0, input.x1), context)

    solution = solve(problem, intctx.solver;
                     merge(input.solver_kwargs,
                          (;callback=CallbackSet(save_callback, update_callback),
                            dense=true))...)

    # Accepted-step saves define the ODE output grid. Populate each SLAB field
    # from those saved states instead of projecting onto the legacy grid.
    _slab_steady_store_saved!(vecs, saved_values, params, ug, tim, bbx, bx, betax)
    stop_t = saved_values.t[end]
    stopped = terminal_state[] !== nothing
    nxtr = length(saved_values.t) + (stopped ? 0 : 1)

    final_state = stopped ? terminal_state[] : context.base
    return (base=final_state, x=stop_t, dx=last_dt[],
        nxtr=nxtr, reference=reference[], controls=controls[],
        stopped=stopped, ode_solution=solution, saved_values=saved_values)
end

"""Copy accepted ODE saves into the variable-length SLAB vector container."""
function _slab_steady_store_saved!(vecs::SLAB_Vecs, saved_values, params,
        ug, tim, bbx, bx, betax)
    states =saved_values.saveval
    n = length(saved_values.t)
    for field in fieldnames(typeof(vecs))
        resize!(getfield(vecs, field), n)
    end
    vecs.x .= saved_values.t
    vecs.xccp .= saved_values.t
    vecs.zc .= getproperty.(states, :zc)
    vecs.h .= getproperty.(states, :h)
    vecs.bb .= getproperty.(states, :bb)
    vecs.b .= getproperty.(states, :b)
    vecs.bbx .= bbx
    vecs.bx .= bx
    vecs.cv .= getproperty.(states, :cv)
    vecs.rho .= getproperty.(states, :rho)
    vecs.t .= getproperty.(states, :t)
    vecs.u .= getproperty.(states, :u)
    vecs.uab .= getproperty.(states, :uab)
    vecs.cm .= getproperty.(states, :cm)
    vecs.cmev .= getproperty.(states, :cmev)
    vecs.cmda .= (1 .- vecs.cm) .* params.met.cmdaa
    vecs.cmw .= getproperty.(states, :cmw)
    vecs.cmwv .= getproperty.(states, :cmwv)
    vecs.wc .= getproperty.(states, :wc)
    vecs.vg .= getproperty.(states, :vg)
    vecs.ug .= ug
    vecs.w .= getproperty.(states, :w)
    vecs.v .= getproperty.(states, :v)
    vecs.vx .= getproperty.(states, :vx)
    vecs.tim .= tim
    vecs.beta .= getproperty.(states, :beta)
    vecs.qint .= getproperty.(states, :qint)
    vecs.betax .= betax
    vecs.tccp .= (2 .* vecs.qint) ./ params.spl.qs
    return vecs
end

"""
Use the ODE solution's dense interpolation for the integrated cloud-center
height; all other fields use the same Akima interpolations as the legacy path.
"""
function _slab_output_interpolations(cc::SLAB_CC_Vecs,
        steady_solution::SLAB_ODE_Steady_Solution{I,F,A,O,V}, params) where {I,F,A,O,V}
    legacy = _slab_legacy_interpolations(cc)
    solution = steady_solution.ode_solution
    saved = steady_solution.saved_values
    x = saved.t
    states = saved.saveval

    # Recreate the legacy concentration post-processing at the callback's
    # accepted states. Width and cloud-duration fields remain shared SLAB
    # outputs; the thermodynamic inputs come directly from the saved ODE states.
    bbx = legacy.bbx_x.(x)
    bx = legacy.bx_x.(x)
    tcld = legacy.tcld.(x)
    points = map(eachindex(x)) do i
        state = states[i]
        _slab_editcc_point(params,x[i],x[1],state.zc,state.h,state.b,state.beta,
            state.uab,state.cm,state.cv,bx[i],bbx[i],tcld[i])
    end

    # Append any transient samples after the callback saves, preserving one
    # simple Akima interpolation per field over the complete output domain.
    tail = findall(xq -> xq > last(x), cc.x)
    x = vcat(x, cc.x[tail])
    centerline = vcat(getproperty.(points, :cc), cc.cc[tail])
    width = vcat(getproperty.(states, :b), cc.b[tail])
    meandered_beta = vcat(getproperty.(points, :betac), cc.betac[tail])
    height = vcat(getproperty.(states, :zc), cc.zc[tail])
    dispersion = vcat(getproperty.(points, :sig), cc.sig[tail])
    return SLAB_Interpolations(
        AkimaInterpolation(centerline,x), AkimaInterpolation(width,x),
        AkimaInterpolation(meandered_beta,x), AkimaInterpolation(height,x),
        AkimaInterpolation(dispersion,x), legacy.xc, legacy.bx,
        legacy.betax, legacy.bx_x, legacy.bbx_x, legacy.tcld, solution)
end