"""
Legacy fixed-step RK4 backend for the steady-state phase.

The legacy backend retains its fixed spatial grid and nested substep loop. It
has no continuous trajectory to retain, so its output interpolation is built
from the stored SLAB vectors using Akima interpolation.
"""

"""
Initialize the legacy backend.

The original RK4 algorithm has no persistent integrator state, so initialization
only returns the dispatch token used by the subsequent step calls.
"""
function _slab_steady_integrator(::SLABLegacySolver, I)
    return SLABLegacyIntegrator()
end

"""
Run the original nested fixed-grid steady-state iteration.

The three inner RK4 steps are legacy behavior: each one uses the preceding
step's refreshed reference values and controls before the next output point is
stored.
"""
function _slab_steady_integrate!(integrator::SLABLegacyIntegrator, vecs,
        vars, params, input::SLAB_Steady_IntegratorInput;
        ug, tim, bbx, bx, betax)
    base = input.base
    reference = _slab_steady_reference_state(base)
    controls = input.controls
    x = input.x0
    nxtr = vars.nxi
    dx = (vars.gam - 1) *
        (params.fld.xffm - vecs.x[vars.msfm]) /
        (vars.gam^(params.xtra.nssm*vars.mnfm) - 1)
    nssm = params.xtra.nssm
    cmdaa = params.met.cmdaa
    qs = params.spl.qs
    rhoa = params.met.rhoa
    work = SLAB_Steady_Workspace(zeros(eltype(vecs.x),11),
        zeros(eltype(vecs.x),11), zeros(eltype(vecs.x),11),
        zeros(eltype(vecs.x),3), zeros(eltype(vecs.x),4))
    stopped = false

    # Preserve the legacy geometric grid: each stored point follows nssm RK4
    # substeps, with the phase reference and controls refreshed after every one.
    for nx in vars.nxi:vars.mffm
        for _ in 1:nssm
            xn = x + dx
            step_input = SLAB_Steady_IntegratorInput(params, base, input.idpf,
                x, xn, controls, input.solver_kwargs)
            step_state = SLAB_Steady_StepState(reference, ug)
            step = SLAB_Steady_StepInput(step_input, step_state, work)
            result = _slab_steady_step!(integrator, step)
            base = _slab_steady_loop_state(result)
            x = xn
            reference, controls = _slab_steady_reference_update(result, controls, rhoa)
            dx *= vars.gam
        end

        # Store one output row after its group of inner RK4 steps, and stop the
        # steady phase at the same heat-release threshold as the original code.
        _slab_sub_store!(vecs,nx,x,base.bb,base.b,base.vg,base.cm,base.t,base.rho,
            base.u,base.h,base.cv,base.beta,base.w,base.v,cmdaa,base.cmw,base.cmwv,
            base.cmev,base.uab,base.wc,base.zc,base.qint,tim,bbx,bx,betax,ug,base.vx)
        vecs.tccp[nx] = (base.qint+base.qint)/qs

        if base.qint >= 0.5*params.spl.qtcs
            nxtr = nx
            stopped = true
            break
        end
    end

    return (base=base, x=x, dx=dx, nxtr=nxtr, reference=reference,
        controls=controls, stopped=stopped, ode_solution=nothing, saved_values=nothing)
end

"""
Build the legacy output interpolations from tabulated cloud data.

Spatial fields are sorted and interpolated against `x`; cloud-center and
crosswind fields are sorted and interpolated against time.
"""
function _slab_legacy_interpolations(cc::SLAB_CC_Vecs)
    xperm = sortperm(cc.x)
    tperm = sortperm(cc.t)
    return SLAB_Interpolations(
        AkimaInterpolation(cc.cc[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.b[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.betac[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.zc[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.sig[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.xc[tperm], cc.t[tperm]),
        AkimaInterpolation(cc.bx[tperm], cc.t[tperm]),
        AkimaInterpolation(cc.betax[tperm], cc.t[tperm]),
        AkimaInterpolation(cc.bx[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.bbx[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.tcld[xperm], cc.x[xperm]),nothing)
end

_slab_output_interpolations(cc::SLAB_CC_Vecs, steady_solution, params) =
    _slab_legacy_interpolations(cc)

"""Construct the backwards-compatible output wrapper using legacy interpolations."""
SLAB_Output(params, state, cc) = SLAB_Output(params, state, cc,
    _slab_legacy_interpolations(cc), nothing, nothing)

"""
Advance one legacy steady-state substep with classical four-stage RK4.

The substep evaluates the legacy slope, solve, thermodynamic, evaluation, and
entrainment routines at each RK stage. Its scratch arrays are reused through
`step.workspace`; the final full phase state and derived cloud volume are
returned together.
"""
function _slab_steady_step!(::SLABLegacyIntegrator, step::SLAB_Steady_StepInput)
    input = step.input
    base, params, idpf = input.base, input.params, input.idpf
    x, xf = input.x0, input.x1
    controls = input.controls
    state = step.state.reference
    work = step.workspace
    f = work.f
    sums = work.sum
    dy = work.dy
    dxxi = work.dxxi
    dxrk = work.dxrk
    fill!(sums, zero(eltype(sums)))
    dx = xf - x
    dxxi .= (dx / 2, dx / 2, dx)
    dxrk .= (dx / 6, dx / 3, dx / 3, dx / 6)

    bbv = bv = qint = zc = r = g = gw = sft = sfu = sfy = sfz = zero(eltype(f))
    cm = base.cm
    cv = cmw = cmwv = cmev = t = _cp = zero(eltype(f))
    rho = base.rho
    u, uab, b, bb, beta, h = base.u, base.uab, base.b, base.bb, base.beta, base.h
    vg, vg0, wc, htp = base.vg, state.vg, base.wc, state.htp
    w, v, vx = base.w, base.v, base.vx
    ubs2, fug = zero(eltype(f)), zero(eltype(f))
    ft, fu, fv, fw = base.ft, base.fu, base.fv, base.fw
    xn = x

    # Classical RK4: accumulate the weighted slope at each stage, then pass
    # the resulting increment through SLAB's solve/thermo/eval/entrainment
    # pipeline to prepare state for the next stage.
    for k in 1:4
        _slab_sub_slope!(f, params, rho, x, h, v, w, b,
                 bb, vg, u, wc, cm, ft, fu, fv, fw,
                 params.othr.bse)
        for j in 1:11
            sums[j] += dxrk[k] * f[j]
            if k == 4
                dy[j] = sums[j]
            else
                dy[j] = dxxi[k] * f[j]
                xn = x + dxxi[k]
            end
        end
        bbv,bv,qint,zc,r,g,gw,sft,sfu,sfy,sfz = _slab_sub_solve(
            params, dy, state.bbv, state.bv, state.zc, state.r, h, state.qint)
        cm,cv,cmw,cmwv,cmev,t,rho,_cp = _slab_sub_thermo(
            params, idpf, xn, zero(eltype(f)), controls.rmi, state.t, state.cmev,
            state.cm, state.cmw, state.cmwv, state.cp, r, state.r, sft,
            params.othr.bse)
        u,uab,b,bb,beta,h,zc,vg,vg0,wc,htp = _slab_sub_eval(
            params, xn, controls.alfg, controls.sru0, zc, state.h, state.u, state.uab,
            state.b, state.bb, r, state.r, bv, state.bv, bbv, state.bbv,
            rho, state.rho, vg0, state.wc, cm, state.htp, b, htp,
            uab, state.beta, vg, wc, h, u, bb, sfu, sfz, sfy, g, gw, params.othr.bse)
        w,v,vx,ubs2,fug,ft,fu,fv,fw = _slab_sub_entran(
            params, idpf, xn, zero(eltype(f)), zero(eltype(f)), zero(eltype(f)),
            state.ubs2, u, step.state.ug, vg, uab, rho, zc, t, h, htp, bb, controls.bbx,
            wc, _cp, params.othr.tgon, params.othr.bse, params.othr.urf,
            params.othr.rcf, params.othr.afa)
    end

    # Return the final full phase state separately from cloud volume, which is
    # derived during projection but stored alongside the state by the caller.
    next = SLAB_Steady_Phase_State(r,bbv,bv,g,gw,sft,sfu,sfy,sfz,zc,qint,h,b,bb,rho,t,
        u,uab,beta,vg0,vg,wc,htp,w,v,vx,cm,cmw,cmwv,cmev,_cp,cv,ft,fu,fv,fw,fug,ubs2)
    return next
end
