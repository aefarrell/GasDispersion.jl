"""Persistent OrdinaryDiffEq driver for the steady-state phase."""

function _slab_int_steady_state!(vecs, vars, params, idpf, nxtr,
                                 ::SLABLegacySolver; solver_kwargs=(;))
    return _slab_int_steady_state_legacy!(vecs, vars, params, idpf, nxtr;
                                          solver_kwargs=solver_kwargs)
end

function _slab_int_steady_state!(vecs, vars, params, idpf, nxtr, solver;
                                 solver_kwargs=(;))
    return _slab_int_steady_state_ode!(vecs, vars, params, idpf, nxtr, solver;
                                       solver_kwargs=solver_kwargs)
end

function _slab_int_steady_state_legacy!(vecs, vars, params, idpf, nxtr;
                                        solver_kwargs=(;))
    return _slab_int_steady_state_impl!(vecs, vars, params, idpf, nxtr;
                                         solver=SLABLegacySolver(), solver_kwargs=solver_kwargs)
end

function _slab_int_steady_state_ode!(vecs, vars, params, idpf, nxtr, solver;
                                     solver_kwargs=(;))
    return _slab_int_steady_state_impl!(vecs, vars, params, idpf, nxtr;
                                          solver=solver,
                                          solver_kwargs=solver_kwargs)
end

function _slab_steady_ode_integrator(base::SLAB_Steady_Phase_State, params,
                                     idpf, x, xf, alg; rmi, alfg, sru0, bbx,
                                     solver_kwargs=(;))
    context = _slab_steady_context(params, base, idpf, nothing;
        rmi=rmi, alfg=alfg, sru0=sru0, bbx=bbx, x0=x)
    problem = ODEProblem(_slab_steady_rhs, context.y0, (x, xf), context)
    integrator = init(problem, alg; solver_kwargs...)
    return integrator, context
end

function _slab_steady_integrator(solver, base, params, idpf, x, xf; kwargs...)
    return _slab_steady_ode_integrator(base, params, idpf, x, xf, solver; kwargs...)
end

function _slab_steady_ode_step!(integrator, context, base, x, xf;
                                rmi, alfg, sru0, bbx, solver_kwargs=(;))
    context.base = base
    context.x0 = x
    context.rmi = rmi
    context.alfg = alfg
    context.sru0 = sru0
    context.bbx = bbx
    context.y0 = SVector{11,typeof(base.r)}(
        base.r, base.bbv, base.bv, base.g, base.sft, base.sfu,
        base.gw, base.zc, base.qint, base.sfy, base.sfz)
    reinit!(integrator, context.y0; t0=x, tf=xf, erase_sol=false)
    requested_dt = get(solver_kwargs, :dt, xf - x)
    set_proposed_dt!(integrator, min(requested_dt, xf - x))
    while integrator.t < xf
        step!(integrator)
    end
    return _slab_steady_project(integrator.u, context, xf)
end

function _slab_steady_step!(integrator_context::Tuple, base, params, idpf, x, xf;
                            solver_kwargs=(;), work, state, rmi, alfg, sru0, bbx)
    integrator, context = integrator_context
    return _slab_steady_ode_step!(integrator, context, base, x, xf;
                                  rmi=rmi, alfg=alfg, sru0=sru0, bbx=bbx,
                                  solver_kwargs=solver_kwargs)
end

function _slab_steady_ode_segment(integrator_context::Tuple, x0, x1)
    _, context = integrator_context
    return SLAB_Steady_ODESegment(x0, x1, deepcopy(context))
end